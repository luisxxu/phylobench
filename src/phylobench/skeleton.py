"""Edge distance on the G-only skeleton of a backbone.

The skeleton is the backbone induced on a keep set (the WoL genome leaves). A placement
on a non-G subtree counts on the skeleton edge that subtree hangs from; that is exact
except at polytomies and above the G root, where it counts at the root vertex. Truth is
computed on the skeleton directly, so it needs no map.
"""

import warnings

import numpy as np
import treeswift

from phylobench.accuracy import (
    Backbone,
    canonical,
    contract_bifurcating_root,
    edge_distance,
    placed_nodes,
    signatures,
)

MIN_LEAVES = 4


def skeleton_newick(newick: str, keep: set[str]) -> str:
    """The tree induced on `keep`; unary nodes are suppressed and lengths summed."""
    tree = treeswift.read_tree_newick(newick)
    with warnings.catch_warnings():
        # Internal labels are support values or taxa; the skeleton drops them by design.
        warnings.simplefilter("ignore", UserWarning)
        skeleton = tree.extract_tree_with(keep)
    if sum(1 for _ in skeleton.traverse_leaves()) < MIN_LEAVES:
        raise ValueError(f"skeleton needs at least {MIN_LEAVES} leaves")
    return skeleton.newick()


def skeleton_map(newick: str, skeleton: Backbone) -> np.ndarray:
    """Skeleton node per backbone node, both in `prepare_backbone` preorder.

    A node with a non-trivial G split takes that split's skeleton node. A node with no
    G below it (non-G pendant) or all of G below it (above the G root) takes its
    parent's, so the root is node 0.
    """
    tree = treeswift.read_tree_newick(newick)
    contract_bifurcating_root(tree.root)
    nodes = list(tree.traverse_preorder())
    index = {node: i for i, node in enumerate(nodes)}
    sig = signatures(tree, skeleton.hashes)
    node_map = np.zeros(len(nodes), dtype=np.int64)
    for i, node in enumerate(nodes[1:], start=1):
        key = canonical(sig[node], skeleton.total)
        node_map[i] = skeleton.key_to_node[key] if key else node_map[index[node.parent]]
    return node_map


def skeleton_scores(
    doc: dict,
    backbone: Backbone,
    skeleton: Backbone,
    node_map: np.ndarray,
    true_node: np.ndarray,
    true_weight: np.ndarray,
    queries: list[str],
) -> list[dict]:
    """One {query, edge_dist} per query; unplaced queries are NaN."""
    rows, placed = placed_nodes(doc, backbone, queries)
    edge_dist, _, _ = edge_distance(
        skeleton, node_map[placed], true_node[rows], true_weight[rows]
    )
    scores = np.full(len(queries), np.nan)
    scores[rows] = edge_dist
    return [
        {"query": q, "edge_dist": d}
        for q, d in zip(queries, scores.tolist(), strict=True)
    ]
