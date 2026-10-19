"""Turn DIPPER's backbone-plus-queries Newick into a jplace doc for `score_jplace`.

The doc is NaN-bearing (lengths are unknown): written to disk only for evaluation.
Never feed it to guppy or any strict-JSON consumer.

Edges are matched by canonical splits, so DIPPER's rooting is irrelevant. Caveats:
- Joint placement is documented, not fixed: a query-only clade maps to the backbone
  edge where the clade attaches.
- A query attached at a backbone vertex is placed on the edge above that vertex (the
  first child edge for the root), so it can score up to 1 against a child-edge truth.
"""

import math

import treeswift

from phylobench.accuracy import jplace_edge_nodes, prepare_backbone, true_positions
from phylobench.jplace import numbered_doc


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
    nodes, _ = true_positions(dipper_newick, present, backbone)

    doc = numbered_doc(backbone_newick)
    edge_of: dict[int, int] = {}
    for number, node in jplace_edge_nodes(doc["tree"], backbone).items():
        edge_of[node] = min(number, edge_of.get(node, number))
    # Root vertex has no edge of its own; edge 1 is the first child edge of the root.
    edge_of[0] = edge_of[1]

    doc["placements"] = [
        {"p": [[edge_of[node], 1.0, math.nan, math.nan]], "n": [query]}
        for query, node in zip(present, nodes.tolist(), strict=True)
    ]
    return doc
