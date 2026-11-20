import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tree_helpers import place

from phylobench.jplace import edge_numbered_tree, numbered_doc

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "skeleton_table.py"
TINY = ROOT / "tests" / "fixtures" / "tiny"
CONFIG = """\
datasets: [gg2-wol]
backbones: [Bfull]
regions: [full-length]
alns: [hmmalign]
exclude:
  - {arm: [pplacer, epa-ng, dipper-r, dipper-m, dipper-m-gpu, apples2, sepp-pinned,
           pplacer-scampp, epa-ng-scampp, bscampp-b5000]}
"""
ACCURACY_HEADER = "run_id\tquery\tsim\tbin\tedge_dist\texact\tdelta_err\tlabel\n"


def genome_names(text: str) -> str:
    """t6.. become G leaves; t1..t5 stay non-G."""
    return re.sub(
        r"\bt(\d+)\b", lambda m: f"G{int(m[1]):09d}" if int(m[1]) > 5 else m[0], text
    )


def build(root: Path):
    data = root / "data" / "gg2-wol" / "Bfull"
    data.mkdir(parents=True)
    backbone = genome_names((TINY / "tree.nwk").read_text())
    (data / "tree.nwk").write_text(backbone)
    (data / "reference_tree.nwk").write_text(
        genome_names((TINY / "reference_tree.nwk").read_text())
    )
    queries = [
        genome_names(line[1:].strip())
        for line in (TINY / "queries.fasta").read_text().splitlines()
        if line.startswith(">")
    ]
    regions = root / "results" / "gg2-wol" / "regions" / "full-length"
    regions.mkdir(parents=True)
    (regions / "queries.fasta").write_text("".join(f">{q}\nACGT\n" for q in queries))
    doc = numbered_doc(backbone)
    n_edges = len(edge_numbered_tree(doc["tree"])[1])
    for i, query in enumerate(queries):
        place(doc, query, i * 7 % n_edges)
    return doc, queries


def write_cell(root: Path, aln: str, arm: str, queries, doc=None):
    cell = root / "results" / "gg2-wol" / "Bfull" / "full-length" / aln / arm
    cell.mkdir(parents=True)
    run_id = f"gg2-wol/Bfull/full-length/{aln}/{arm}"
    rows = "".join(
        f"{run_id}\t{q}\t0.9\t90-95\t0.0\t1.0\t0.0\taccuracy\n" for q in queries
    )
    (cell / "accuracy.tsv").write_text(ACCURACY_HEADER + (rows if doc else ""))
    if doc:
        (cell / "placements.jplace").write_text(json.dumps(doc))


def run(root: Path):
    (root / "t.yaml").write_text(CONFIG)
    return subprocess.run(
        [sys.executable, str(SCRIPT), "t.yaml"],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        cwd=root,
    )


def test_scores_ok_cells_including_sepp_prod_and_skips_dnf_cells(tmp_path):
    doc, queries = build(tmp_path)
    write_cell(tmp_path, "raw", "sepp-prod", queries, doc)
    write_cell(tmp_path, "hmmalign", "bscampp", queries)
    result = run(tmp_path)
    assert result.returncode == 0, result.stderr
    table = pd.read_csv(tmp_path / "results" / "skeleton.tsv", sep="\t")
    assert list(table.columns) == ["run_id", "query", "edge_dist"]
    assert set(table["run_id"]) == {"gg2-wol/Bfull/full-length/raw/sepp-prod"}
    assert sorted(table["query"]) == sorted(queries)
    assert table["edge_dist"].notna().all()
    assert "1 scored, 1 skipped" in result.stdout
