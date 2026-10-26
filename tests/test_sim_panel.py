import runpy
import sys
from pathlib import Path

import pytest
import treeswift
from emp_helpers import run_script

from phylobench.accuracy import prepare_backbone, true_positions

SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "make_sim_panel.py"
N_QUERIES = 7 * (5 + 10 + 15 + 20 + 25) + 500


def read_fasta_names(path):
    return [line[1:] for line in path.read_text().splitlines() if line[0] == ">"]


@pytest.fixture(scope="module")
def panel(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("panel")
    argv, sys.argv = sys.argv, ["make_sim_panel.py", str(out_dir)]
    try:
        runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
    return out_dir


def test_sim_panel_meets_acceptance_properties(panel):
    tmp_path = panel

    reference_text = (tmp_path / "reference_tree.nwk").read_text()
    backbone_text = (tmp_path / "tree.nwk").read_text()
    truth = [
        line.split("\t")
        for line in (tmp_path / "truth.tsv").read_text().splitlines()[1:]
    ]
    queries = [row[0] for row in truth]
    assert len(queries) == len(set(queries)) == N_QUERIES
    assert read_fasta_names(tmp_path / "queries.aligned.fasta") == queries

    backbone = treeswift.read_tree_newick(backbone_text)
    leaves = {leaf.label for leaf in backbone.traverse_leaves()}
    assert len(leaves) == 10_000 - N_QUERIES
    assert leaves.isdisjoint(queries)
    assert set(read_fasta_names(tmp_path / "aligned-dna-sequences.fasta")) == leaves
    assert all(len(n.children) != 1 for n in backbone.traverse_preorder())
    assert len(backbone.root.children) != 2

    nodes, _ = true_positions(reference_text, queries, prepare_backbone(backbone_text))
    assert len(nodes) == N_QUERIES


def test_every_sim_query_survives_v4_150_extraction(panel, tmp_path):
    result = run_script(
        "extract_regions.py",
        "--region", "V4-150",
        "--queries", panel / "queries.fasta",
        "--out-fasta", tmp_path / "queries.fasta",
        "--out-exclusions", tmp_path / "exclusions.tsv",
        "--out-counts", tmp_path / "counts.tsv",
    )  # fmt: skip
    assert result.returncode == 0, result.stderr
    assert read_fasta_names(tmp_path / "queries.fasta") == read_fasta_names(
        panel / "queries.fasta"
    )
    assert (tmp_path / "exclusions.tsv").read_text() == ""
