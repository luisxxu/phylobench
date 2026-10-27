"""Turn DIPPER's backbone-plus-queries Newick into a jplace doc for `score_jplace`.

The doc is NaN-bearing (lengths are unknown): written to disk only for evaluation.
Never feed it to guppy or any strict-JSON consumer.

Edges are matched by canonical splits, so DIPPER's rooting is irrelevant. Caveats:
- Joint placement is documented, not fixed: a query-only clade maps to the backbone
  edge where the clade attaches.
- A query attached at a backbone vertex is placed on the edge above that vertex (the
  first child edge for the root), so it can score up to 1 against a child-edge truth.
- Zero-length backbone internal edges are contracted for matching (CPU dipper-m always
  collapses them); a query at such a merged vertex lands on the edge above it.
"""

import math

import treeswift

from phylobench.accuracy import (
    ZERO_LENGTH,
    Backbone,
    canonical,
    contract_bifurcating_root,
    jplace_edge_nodes,
    prepare_backbone,
    signatures,
    true_positions,
)
from phylobench.jplace import numbered_doc


def _zero_length_splits(backbone_newick: str, backbone: Backbone) -> set[int]:
    tree = treeswift.read_tree_newick(backbone_newick)
    contract_bifurcating_root(tree.root)
    sig = signatures(tree, backbone.hashes)
    return {
        canonical(sig[node], backbone.total)
        for node in tree.traverse_internal()
        if node is not tree.root
        and node.edge_length is not None
        and node.edge_length <= ZERO_LENGTH
    }


def _contract(newick: str, backbone: Backbone, splits: set[int]) -> str:
    """Newick with each internal edge whose backbone split is in `splits` contracted."""
    tree = treeswift.read_tree_newick(newick)
    sig = signatures(tree, backbone.hashes)
    for node in list(tree.traverse_preorder()):
        if (
            node is not tree.root
            and not node.is_leaf()
            and canonical(sig[node], backbone.total) in splits
        ):
            siblings = node.parent.children
            at = siblings.index(node)
            siblings[at : at + 1] = node.children
            for child in node.children:
                child.parent = node.parent
    return tree.newick()


def dipper_to_jplace(
    dipper_newick: str, backbone_newick: str, queries: list[str]
) -> dict:
    """jplace doc with each query on the backbone edge DIPPER attached it to.

    Queries absent from DIPPER's tree get no record (NaN row when scored). Raises on
    duplicate queries, a query named like a backbone leaf, or any difference between
    DIPPER's tree restricted to the backbone and the backbone.
    """
    if len(set(queries)) != len(queries):
        raise ValueError("duplicate entries in queries")
    backbone = prepare_backbone(backbone_newick)
    clashes = backbone.hashes.keys() & set(queries)
    if clashes:
        raise ValueError(f"queries share names with backbone leaves: {sorted(clashes)}")

    dipper_leaves = {
        leaf.label
        for leaf in treeswift.read_tree_newick(dipper_newick).traverse_leaves()
    }
    present = [q for q in queries if q in dipper_leaves]
    if dipper_leaves != backbone.hashes.keys() | set(present):
        raise ValueError("DIPPER leaves are not exactly the backbone plus the queries")
    zero_splits = _zero_length_splits(backbone_newick, backbone)
    if zero_splits:
        contracted = prepare_backbone(_contract(backbone_newick, backbone, zero_splits))
        contracted_key = {i: k for k, i in contracted.key_to_node.items()}
        contracted_nodes, _ = true_positions(
            _contract(dipper_newick, backbone, zero_splits), present, contracted
        )
        nodes = [
            backbone.key_to_node[contracted_key[i]] if i else 0
            for i in contracted_nodes.tolist()
        ]
    else:
        nodes = true_positions(dipper_newick, present, backbone)[0].tolist()

    doc = numbered_doc(backbone_newick)
    edge_of: dict[int, int] = {}
    for number, node in jplace_edge_nodes(doc["tree"], backbone).items():
        edge_of[node] = min(number, edge_of.get(node, number))
    # Root vertex has no edge of its own; edge 1 is the first child edge of the root.
    edge_of[0] = edge_of[1]

    doc["placements"] = [
        {"p": [[edge_of[node], 1.0, math.nan, math.nan]], "n": [query]}
        for query, node in zip(present, nodes, strict=True)
    ]
    return doc
