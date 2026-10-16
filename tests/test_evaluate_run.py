import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from tree_helpers import place

from phylobench.accuracy import jplace_edge_nodes, prepare_backbone, true_positions
from phylobench.jplace import numbered_doc
from phylobench.regions import read_fasta
from phylobench.tables import ACCURACY_COLUMNS

ROOT = Path(__file__).parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
FIXTURE = ROOT / "tests" / "fixtures" / "tiny"
QUERIES = [name for name, _ in read_fasta(FIXTURE / "queries.fasta")]


def perfect_doc():
    tree = (FIXTURE / "tree.nwk").read_text()
    backbone = prepare_backbone(tree)
    doc = numbered_doc(tree)
    edge_of = {
        node: edge for edge, node in jplace_edge_nodes(doc["tree"], backbone).items()
    }
    nodes, _ = true_positions(
        (FIXTURE / "reference_tree.nwk").read_text(), QUERIES, backbone
    )
    for query, node in zip(QUERIES, nodes.tolist(), strict=True):
        place(doc, query, edge_of[node])
    return doc


def run_script(script, *args):
    return subprocess.run(
        [sys.executable, str(SCRIPTS / script), *map(str, args)],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )


@pytest.fixture
def cell(tmp_path):
    (tmp_path / "region.fasta").write_text("".join(f">{q}\nACGT\n" for q in QUERIES))
    pd.DataFrame({"query": QUERIES, "sim": 0.9, "bin": "90-95"}).to_csv(
        tmp_path / "sim.tsv", sep="\t", index=False
    )
    return tmp_path


def evaluate(cell, status, jplace_text):
    (cell / "gate.run.json").write_text(json.dumps({"status": status}))
    (cell / "placements.jplace").write_text(jplace_text)
    return run_script(
        "evaluate_run.py",
        "--run-id", "rid",
        "--gate-run-json", cell / "gate.run.json",
        "--jplace", cell / "placements.jplace",
        "--backbone-tree", FIXTURE / "tree.nwk",
        "--reference-tree", FIXTURE / "reference_tree.nwk",
        "--region-queries", cell / "region.fasta",
        "--similarity", cell / "sim.tsv",
        "--out", cell / "accuracy.tsv",
    )  # fmt: skip


def test_ok_run_scores_every_query_exact(cell):
    result = evaluate(cell, "ok", json.dumps(perfect_doc()))
    assert result.returncode == 0, result.stderr
    table = pd.read_csv(cell / "accuracy.tsv", sep="\t")
    assert tuple(table.columns) == ACCURACY_COLUMNS
    assert sorted(table["query"]) == sorted(QUERIES)
    assert table["exact"].sum() == 10
    assert (table["sim"] == 0.9).all() and (table["label"] == "accuracy").all()


def test_timeout_gate_writes_header_only_without_reading_jplace(cell):
    result = evaluate(cell, "timeout", "")
    assert result.returncode == 0, result.stderr
    table = pd.read_csv(cell / "accuracy.tsv", sep="\t")
    assert table.empty and tuple(table.columns) == ACCURACY_COLUMNS


def test_ok_gate_with_empty_jplace_fails(cell):
    result = evaluate(cell, "ok", "")
    assert result.returncode != 0
    assert "empty jplace" in result.stderr


def test_similarity_script_one_row_per_query_in_unit_interval(tmp_path):
    out = tmp_path / "similarity.tsv"
    result = run_script(
        "similarity.py",
        "--queries", FIXTURE / "queries.aligned.fasta",
        "--backbone", FIXTURE / "aligned-dna-sequences.fasta",
        "--out", out,
    )  # fmt: skip
    assert result.returncode == 0, result.stderr
    table = pd.read_csv(out, sep="\t")
    assert list(table.columns) == ["query", "sim", "bin"]
    assert sorted(table["query"]) == sorted(QUERIES)
    assert table["sim"].between(0, 1).all() and (table["sim"] > 0).all()
