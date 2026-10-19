import json
import subprocess
from pathlib import Path

import pandas as pd
from test_evaluate_run import FIXTURE, QUERIES, run_script

from phylobench.regions import read_fasta

ROOT = Path(__file__).parents[1]
UNGAP = ROOT / "workflow" / "scripts" / "ungap_concat.awk"


def to_jplace(tmp_path, status, dipper_tree):
    (tmp_path / "gate.run.json").write_text(json.dumps({"status": status}))
    (tmp_path / "dipper.nwk").write_text(dipper_tree)
    return run_script(
        "dipper_to_jplace.py",
        "--gate-run-json", tmp_path / "gate.run.json",
        "--dipper-tree", tmp_path / "dipper.nwk",
        "--backbone-tree", FIXTURE / "tree.nwk",
        "--region-queries", FIXTURE / "queries.fasta",
        "--out", tmp_path / "placements.jplace",
    )  # fmt: skip


def test_ok_output_scores_every_query_exact(tmp_path):
    result = to_jplace(tmp_path, "ok", (FIXTURE / "reference_tree.nwk").read_text())
    assert result.returncode == 0, result.stderr
    text = (tmp_path / "placements.jplace").read_text()
    assert "NaN" in text
    json.loads(text)

    (tmp_path / "sim.tsv").write_text(
        "query\tsim\tbin\n" + "".join(f"{q}\t0.9\t90-95\n" for q in QUERIES)
    )
    scored = run_script(
        "evaluate_run.py",
        "--run-id", "rid",
        "--gate-run-json", tmp_path / "gate.run.json",
        "--jplace", tmp_path / "placements.jplace",
        "--backbone-tree", FIXTURE / "tree.nwk",
        "--reference-tree", FIXTURE / "reference_tree.nwk",
        "--region-queries", FIXTURE / "queries.fasta",
        "--similarity", tmp_path / "sim.tsv",
        "--out", tmp_path / "accuracy.tsv",
    )  # fmt: skip
    assert scored.returncode == 0, scored.stderr
    table = pd.read_csv(tmp_path / "accuracy.tsv", sep="\t")
    assert sorted(table["query"]) == sorted(QUERIES)
    assert (table["edge_dist"] == 0).all()


def test_dnf_gate_with_empty_tree_writes_empty_jplace(tmp_path):
    result = to_jplace(tmp_path, "timeout", "")
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "placements.jplace").read_text() == ""


def test_changed_backbone_topology_fails(tmp_path):
    tree = (FIXTURE / "reference_tree.nwk").read_text()
    leaves = [name for name, _ in read_fasta(FIXTURE / "aligned-dna-sequences.fasta")]
    a, b = leaves[:2]
    swapped = tree.replace(a, "@@").replace(b, a).replace("@@", b)
    result = to_jplace(tmp_path, "ok", swapped)
    assert result.returncode != 0 and "ValueError" in result.stderr


def awk(*paths):
    return subprocess.run(
        ["awk", "-f", str(UNGAP), *map(str, paths)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout


def test_ungap_concat_strips_gaps_and_crlf_across_files(tmp_path):
    backbone = tmp_path / "backbone.fasta"
    backbone.write_bytes(b">a\r\nAC-GT\r\n..AC\r\n\r\n>b\nA-C.G")
    queries = tmp_path / "queries.fasta"
    queries.write_text(">q1\nACGT\n>q2\nGG\n")
    combined = tmp_path / "combined.fasta"
    combined.write_text(awk(backbone) + awk(queries))

    assert "\r" not in combined.read_text()
    assert "" not in combined.read_text().splitlines()
    assert list(read_fasta(combined)) == [
        ("a", "ACGTAC"),
        ("b", "ACG"),
        ("q1", "ACGT"),
        ("q2", "GG"),
    ]


def test_awk_1_terminates_a_missing_final_newline_and_keeps_blank_lines(tmp_path):
    aligned = tmp_path / "aligned.fasta"
    aligned.write_text(">a\nAC\n\n>b\nGT")
    out = subprocess.run(
        ["awk", "1", str(aligned), str(aligned)], capture_output=True, text=True
    ).stdout
    assert out == ">a\nAC\n\n>b\nGT\n>a\nAC\n\n>b\nGT\n"
