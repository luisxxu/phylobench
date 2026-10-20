"""Polytomy-aware placement accuracy against a true tree.

Edges are identified across trees by canonical XOR-hash splits of backbone leaves, so
topologies need not share a rooting. Zero-length backbone edges are collapsed: a point
on such an edge is the vertex it collapses into.
"""

from dataclasses import dataclass

import numpy as np
import treeswift

from phylobench.jplace import edge_numbered_tree, top_placements

ZERO_LENGTH = 0.0


def leaf_hashes(labels, seed: int = 0) -> dict[str, int]:
    """Seeded uint64 per label, independent of the order labels are given in."""
    ordered = sorted(labels)
    values = np.random.default_rng(seed).integers(
        1, np.iinfo(np.uint64).max, size=len(ordered), dtype=np.uint64, endpoint=True
    )
    return dict(zip(ordered, values.tolist(), strict=True))


def signatures(
    tree: treeswift.Tree, hashes: dict[str, int]
) -> dict[treeswift.Node, int]:
    """XOR of the hashes below each node; leaves without a hash contribute 0."""
    sig: dict[treeswift.Node, int] = {}
    for node in tree.traverse_postorder():
        if node.is_leaf():
            sig[node] = hashes.get(node.label, 0)
        else:
            value = 0
            for child in node.children:
                value ^= sig[child]
            sig[node] = value
    return sig


def canonical(sig: int, total: int) -> int:
    """Orientation-free split key; 0 means the split is trivial for the backbone."""
    return min(sig, sig ^ total)


@dataclass(frozen=True, eq=False)
class Backbone:
    """Backbone tree flattened to preorder indices (root is 0)."""

    hashes: dict[str, int]
    total: int
    key_to_node: dict[int, int]
    sig: list[int]
    weight: np.ndarray
    nz_depth: np.ndarray
    depth: np.ndarray
    up: np.ndarray


def _sum_lengths(a, b):
    return None if a is None or b is None else a + b


def contract_bifurcating_root(root: treeswift.Node) -> None:
    """Drop a degree-2 root: its internal child merges in; the other edge sums both."""
    if len(root.children) != 2:
        return
    inner, other = root.children
    if inner.is_leaf():
        inner, other = other, inner
    if inner.is_leaf():
        return
    other.edge_length = _sum_lengths(inner.edge_length, other.edge_length)
    root.children = [other, *inner.children]
    for child in inner.children:
        child.parent = root


def prepare_backbone(newick: str, seed: int = 0) -> Backbone:
    tree = treeswift.read_tree_newick(newick)
    contract_bifurcating_root(tree.root)
    nodes = list(tree.traverse_preorder())
    index = {node: i for i, node in enumerate(nodes)}
    hashes = leaf_hashes(tree.labels(internal=False), seed)
    node_sig = signatures(tree, hashes)
    total = node_sig[tree.root]

    count = len(nodes)
    parent = np.zeros(count, dtype=np.int32)
    weight = np.zeros(count, dtype=np.int8)
    depth = np.zeros(count, dtype=np.int32)
    nz_depth = np.zeros(count, dtype=np.int32)
    key_to_node: dict[int, int] = {}
    for i, node in enumerate(nodes[1:], start=1):
        p = index[node.parent]
        length = node.edge_length
        parent[i] = p
        weight[i] = 0 if length is not None and length <= ZERO_LENGTH else 1
        depth[i] = depth[p] + 1
        nz_depth[i] = nz_depth[p] + weight[i]
        key = canonical(node_sig[node], total)
        if key in key_to_node:
            raise ValueError(
                "backbone has duplicate splits (unary node or repeated leaf)"
            )
        key_to_node[key] = i

    levels = max(1, int(depth.max()).bit_length())
    up = np.empty((levels, count), dtype=np.int32)
    up[0] = parent
    for k in range(1, levels):
        up[k] = up[k - 1][up[k - 1]]
    return Backbone(
        hashes,
        total,
        key_to_node,
        [node_sig[n] for n in nodes],
        weight,
        nz_depth,
        depth,
        up,
    )


def jplace_edge_nodes(newick: str, backbone: Backbone) -> dict[int, int]:
    """Map each jplace edge number to its backbone node (preorder index).

    The two root edges of a bifurcating jplace root are one backbone edge, so they map
    to the same node.
    """
    tree, edges = edge_numbered_tree(newick)
    if set(tree.labels(internal=False)) != backbone.hashes.keys():
        raise ValueError("jplace tree leaves differ from the backbone leaves")
    sig = signatures(tree, backbone.hashes)
    root = tree.root
    root_pair = set(root.children) if len(root.children) == 2 else set()

    mapping: dict[int, int] = {}
    first_with_key: dict[int, treeswift.Node] = {}
    for number, node in edges.items():
        if node is root:
            continue
        key = canonical(sig[node], backbone.total)
        if key not in backbone.key_to_node:
            raise ValueError(f"jplace edge {number} is not a backbone edge")
        first = first_with_key.setdefault(key, node)
        if first is not node and not (node in root_pair and first in root_pair):
            raise ValueError(f"jplace edge {number} duplicates another edge's split")
        mapping[number] = backbone.key_to_node[key]
    return mapping


def true_positions(
    reference_newick: str, queries: list[str], backbone: Backbone
) -> tuple[np.ndarray, np.ndarray]:
    """Where each query attaches in the reference restricted to the backbone.

    Returns (node, true_weight). An edge truth is the backbone node below that edge
    with its weight; a vertex truth (three or more backbone directions meet at the
    attachment) is the vertex with weight 0.
    """
    tree = treeswift.read_tree_newick(reference_newick)
    # Tree.label_to_node("leaves") wrongly includes an unlabelled root.
    leaves = {leaf.label: leaf for leaf in tree.traverse_leaves()}
    backbone_labels = backbone.hashes.keys()
    if not backbone_labels.isdisjoint(queries) or not (
        backbone_labels | set(queries) <= leaves.keys()
    ):
        raise ValueError(
            "reference leaves must include all backbone leaves and the queries, "
            "which must not be backbone leaves"
        )
    sig = signatures(tree, backbone.hashes)
    total = backbone.total
    reference_keys = {canonical(s, total) for s in sig.values()} - {0}
    if reference_keys != backbone.key_to_node.keys():
        raise ValueError(
            "reference restricted to the backbone differs from the backbone"
        )

    nearest_carrier: dict[treeswift.Node, treeswift.Node] = {}
    for node in tree.root.traverse_preorder(include_self=False):
        nearest_carrier[node] = (
            node.parent if sig[node.parent] else nearest_carrier[node.parent]
        )

    nodes, weights = [], []
    for query in queries:
        attach = nearest_carrier[leaves[query]]
        while True:
            carriers = [c for c in attach.children if sig[c]]
            directions = len(carriers) + (sig[attach] != total)
            if directions != 1:
                break
            attach = carriers[0]
        child = carriers[0]
        node = backbone.key_to_node[canonical(sig[child], total)]
        if directions == 2:
            weights.append(backbone.weight[node])
        else:
            if backbone.sig[node] == sig[child]:
                node = int(backbone.up[0][node])
            weights.append(0)
        nodes.append(node)
    return np.array(nodes, dtype=np.int64), np.array(weights, dtype=np.int64)


def _lca(up, depth, u, v):
    deeper_is_v = depth[u] < depth[v]
    u, v = np.where(deeper_is_v, v, u), np.where(deeper_is_v, u, v)
    gap = depth[u] - depth[v]
    for k in range(up.shape[0]):
        u = np.where((gap >> k) & 1 == 1, up[k][u], u)
    for k in reversed(range(up.shape[0])):
        above_u, above_v = up[k][u], up[k][v]
        diverged = above_u != above_v
        u, v = np.where(diverged, above_u, u), np.where(diverged, above_v, v)
    return np.where(u == v, u, up[0][u])


def edge_distance(
    backbone: Backbone,
    placed: np.ndarray,
    true_node: np.ndarray,
    true_weight: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (edge_dist, exact, delta_err) per placement, all vectorized.

    edge_dist is the number of collapsed-tree vertices on the path between the placed
    edge and the true edge (0 on the same edge, 1 for adjacent edges); a vertex truth is
    0 from every edge incident to it. delta_err counts reference bipartitions over
    backbone + query missing from the estimate. No NaN arises here; score_jplace owns
    the missing-placement contract.

    T19 concordance calls this with truth = the sepp-prod node and
    true_weight = weight[node], and uses edge_dist only.
    """
    p = np.asarray(placed, dtype=np.int64)
    t = np.asarray(true_node, dtype=np.int64)
    w_t = np.asarray(true_weight, dtype=np.int64)
    nz = backbone.nz_depth.astype(np.int64)
    w = backbone.weight.astype(np.int64)
    lca = _lca(backbone.up, backbone.depth, p, t)

    above_p, above_t = nz[p] - w[p], nz[t] - w_t
    neither = (above_p - nz[lca]) + (above_t - nz[lca])
    between = np.where(
        lca == p, above_t - nz[p], np.where(lca == t, above_p - nz[t], neither)
    )
    same = p == t
    edge_dist = np.where(same, 0, between + w[p] * w_t)
    delta_err = np.where(same, 0, between + w_t)
    return edge_dist, edge_dist == 0, delta_err


def placed_nodes(
    doc: dict, backbone: Backbone, queries: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    """(row indices into queries, backbone nodes) of the placed queries.

    Unplaced queries are absent. Raises on duplicate requested queries, placements of
    queries that were not requested, and jplace leaves that differ from the backbone.
    """
    if len(set(queries)) != len(queries):
        raise ValueError("duplicate entries in queries")
    placements = top_placements(doc)
    unknown = {p.query for p in placements} - set(queries)
    if unknown:
        raise ValueError(
            f"jplace places queries that were not requested: {sorted(unknown)}"
        )
    edge_nodes = jplace_edge_nodes(doc["tree"], backbone)
    position = {q: i for i, q in enumerate(queries)}
    rows = np.array([position[p.query] for p in placements], dtype=np.int64)
    placed = np.array([edge_nodes[p.edge_num] for p in placements], dtype=np.int64)
    return rows, placed


def score_jplace(
    doc: dict, backbone_newick: str, reference_newick: str, queries: list[str]
) -> list[dict]:
    """One row per requested query: {query, edge_dist, exact (1.0/0.0), delta_err}.

    Unplaced queries are all NaN. Consumers must not skipna: count NaN as missing and
    report the count.
    """
    backbone = prepare_backbone(backbone_newick)
    rows, placed = placed_nodes(doc, backbone, queries)
    true_node, true_weight = true_positions(reference_newick, queries, backbone)
    edge_dist, exact, delta_err = edge_distance(
        backbone, placed, true_node[rows], true_weight[rows]
    )

    scores = np.full((len(queries), 3), np.nan)
    scores[rows] = np.column_stack([edge_dist, exact, delta_err])
    return [
        {"query": q, "edge_dist": s[0], "exact": s[1], "delta_err": s[2]}
        for q, s in zip(queries, scores.tolist(), strict=True)
    ]
