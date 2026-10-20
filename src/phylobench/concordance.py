"""Edge distance between two placements of the same queries on one backbone.

Concordance is not accuracy: neither placement is truth. It reuses the T5 core with
the reference placement standing in for the true position.
"""

import numpy as np

from phylobench.accuracy import edge_distance, placed_nodes, prepare_backbone


def score_concordance(
    doc: dict, reference_doc: dict, backbone_newick: str, queries: list[str]
) -> list[dict]:
    """One row per requested query: {query, edge_dist}.

    NaN when either doc leaves the query unplaced. Consumers must not skipna.
    """
    backbone = prepare_backbone(backbone_newick)
    rows, nodes = placed_nodes(doc, backbone, queries)
    reference_rows, reference_nodes = placed_nodes(reference_doc, backbone, queries)
    both, in_doc, in_reference = np.intersect1d(
        rows, reference_rows, return_indices=True
    )
    reference = reference_nodes[in_reference]
    edge_dist, _, _ = edge_distance(
        backbone, nodes[in_doc], reference, backbone.weight[reference]
    )
    scores = np.full(len(queries), np.nan)
    scores[both] = edge_dist
    return [
        {"query": q, "edge_dist": d}
        for q, d in zip(queries, scores.tolist(), strict=True)
    ]
