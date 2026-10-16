import json
import os
import subprocess
from pathlib import Path

import pytest

from phylobench.accuracy import jplace_edge_nodes, prepare_backbone, true_positions
from phylobench.jplace import edge_numbered_tree, numbered_doc

GUPPY = os.environ.get("GUPPY")
TINY = Path(__file__).parent / "fixtures" / "tiny"

pytestmark = pytest.mark.skipif(not GUPPY, reason="set GUPPY to the guppy binary")


def test_tog_inserts_every_query_on_its_placed_edge(tmp_path):
    backbone = (TINY / "tree.nwk").read_text().strip()
    queries = [
        line[1:].strip()
        for line in (TINY / "queries.fasta").read_text().splitlines()
        if line.startswith(">")
    ]
    doc = numbered_doc(
        backbone,
        fields=[
            "edge_num",
            "likelihood",
            "like_weight_ratio",
            "distal_length",
            "pendant_length",
        ],
    )
    doc["metadata"] = {"invocation": "test_guppy_tog"}
    tree, edges = edge_numbered_tree(doc["tree"])
    placeable = sorted(
        n for n, node in edges.items() if node.edge_length and node.edge_length > 0
    )
    assert len(placeable) >= len(queries)
    placed_on = dict(zip(queries, placeable, strict=False))
    for query, edge in placed_on.items():
        length = edges[edge].edge_length
        doc["placements"].append({"p": [[edge, -1.0, 1.0, 0.5 * length, 0.01]], "n": [query]})
    jplace = tmp_path / "in.jplace"
    jplace.write_text(json.dumps(doc))
    out = tmp_path / "out.tre"

    subprocess.run([GUPPY, "tog", "-o", str(out), str(jplace)], check=True)

    inserted = out.read_text()
    prepared = prepare_backbone(backbone)
    nodes, _ = true_positions(inserted, queries, prepared)
    edge_nodes = jplace_edge_nodes(doc["tree"], prepared)
    for query, node in zip(queries, nodes.tolist(), strict=True):
        assert node == edge_nodes[placed_on[query]]
