import json
import math

import pandas as pd
import pytest
from emp_helpers import ROOT, run_script, write, write_run_json
from tree_helpers import edge_numbers, place

from phylobench.jplace import numbered_doc
from phylobench.tables import CONCORDANCE_COLUMNS

TREE = (ROOT / "tests" / "fixtures" / "tiny" / "tree.nwk").read_text()
QUERIES = ["q1", "q2", "q3"]


def jplace_text(placed):
    doc = numbered_doc(TREE)
    numbers = edge_numbers(doc)
    for query, leaf in placed.items():
        place(doc, query, numbers[leaf])
    return json.dumps(doc)


def concordance(tmp_path, gate="ok", reference_gate="ok", jplace=None, reference=None):
    write_run_json(tmp_path / "gate.json", gate)
    write_run_json(tmp_path / "reference.json", reference_gate)
    write(
        tmp_path / "arm.jplace",
        jplace_text({"q1": "t110", "q2": "t109"}) if jplace is None else jplace,
    )
    write(
        tmp_path / "sepp.jplace",
        jplace_text({"q1": "t110", "q2": "t110", "q3": "t110"})
        if reference is None
        else reference,
    )
    write(tmp_path / "queries.fasta", "".join(f">{q}\nACGT\n" for q in QUERIES))
    write(tmp_path / "tree.nwk", TREE)
    pd.DataFrame({"query": QUERIES, "sim": 0.9, "bin": "90-95"}).to_csv(
        tmp_path / "sim.tsv", sep="\t", index=False
    )
    return run_script(
        "concordance_run.py",
        "--run-id", "emp-10k/B10k/V4-150/hmmalign/pplacer",
        "--reference-run-id", "emp-10k/B10k/V4-150/raw/sepp-prod",
        "--gate-run-json", tmp_path / "gate.json",
        "--reference-gate-run-json", tmp_path / "reference.json",
        "--jplace", tmp_path / "arm.jplace",
        "--reference-jplace", tmp_path / "sepp.jplace",
        "--backbone-tree", tmp_path / "tree.nwk",
        "--region-queries", tmp_path / "queries.fasta",
        "--similarity", tmp_path / "sim.tsv",
        "--out", tmp_path / "concordance.tsv",
    )  # fmt: skip


def read_table(tmp_path):
    return pd.read_csv(tmp_path / "concordance.tsv", sep="\t")


def test_ok_cells_write_labelled_concordance_rows(tmp_path):
    result = concordance(tmp_path)
    assert result.returncode == 0, result.stderr
    table = read_table(tmp_path)
    assert tuple(table.columns) == CONCORDANCE_COLUMNS
    assert table["query"].tolist() == QUERIES
    assert table["edge_dist"][:2].tolist() == [0, 1]
    assert math.isnan(table["edge_dist"][2])
    assert (table["label"] == "concordance").all()
    assert (table["reference_run_id"] == "emp-10k/B10k/V4-150/raw/sepp-prod").all()


@pytest.mark.parametrize(
    ("gate", "reference_gate"), [("timeout", "ok"), ("ok", "oom"), ("error", "error")]
)
def test_either_gate_not_ok_writes_header_only_without_reading_jplace(
    tmp_path, gate, reference_gate
):
    result = concordance(
        tmp_path, gate, reference_gate, jplace="", reference="not json"
    )
    assert result.returncode == 0, result.stderr
    table = read_table(tmp_path)
    assert table.empty and tuple(table.columns) == CONCORDANCE_COLUMNS


@pytest.mark.parametrize("which", ["jplace", "reference"])
def test_ok_gate_with_empty_jplace_fails(tmp_path, which):
    result = concordance(tmp_path, **{which: ""})
    assert result.returncode != 0
    assert "empty jplace" in result.stderr
