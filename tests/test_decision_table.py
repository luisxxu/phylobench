import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "decision_table.py"
QUERIES = [f"G{i:09d}" for i in range(8)]
CONFIG = """\
datasets: [gg2-wol]
backbones: [Bfull]
regions: [full-length, V4-150]
alns: [hmmalign]
exclude:
  - {arm: [pplacer, epa-ng, dipper-r, dipper-m, dipper-m-gpu, pplacer-scampp,
           epa-ng-scampp, bscampp-b5000]}
"""


def cells(arm, aln, edge_dist):
    return pd.concat(
        [
            pd.DataFrame(
                {
                    "run_id": f"gg2-wol/Bfull/{region}/{aln}/{arm}",
                    "query": QUERIES,
                    "edge_dist": edge_dist,
                }
            )
            for region in ("full-length", "V4-150")
        ]
    )


def run(tmp_path):
    good, bad = [0.0] * 8, [9.0] * 8
    accuracy = pd.concat(
        [
            cells("sepp-prod", "raw", good),
            cells("sepp-pinned", "raw", good),
            cells("bscampp", "hmmalign", bad),
            cells("apples2", "hmmalign", good),
        ]
    )
    accuracy.to_csv(tmp_path / "accuracy.tsv", sep="\t", index=False)
    accuracy.assign(edge_dist=accuracy["edge_dist"].clip(upper=1)).to_csv(
        tmp_path / "skeleton.tsv", sep="\t", index=False
    )
    pd.DataFrame({"leaf": QUERIES, "group": QUERIES, "query": 1}).to_csv(
        tmp_path / "heldout.tsv", sep="\t", index=False
    )
    (tmp_path / "t30.yaml").write_text(CONFIG)
    return subprocess.run(
        [
            sys.executable, str(SCRIPT), str(tmp_path / "t30.yaml"),
            "--accuracy", str(tmp_path / "accuracy.tsv"),
            "--skeleton", str(tmp_path / "skeleton.tsv"),
            "--heldout", str(tmp_path / "heldout.tsv"),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )  # fmt: skip


def test_headline_per_core_arm_and_none_for_cost(tmp_path):
    result = run(tmp_path)
    assert result.returncode == 0, result.stderr
    decision, headline = result.stdout.split("headline", 1)
    assert "COST (descriptive)" in decision and "descriptive" in decision
    lines = [line.split() for line in headline.strip().splitlines()]
    assert {"sepp-pinned", "bscampp"} <= {w for line in lines for w in line}
    assert "apples2" not in headline
    assert "as accurate as production SEPP" in headline
    assert "not shown" in headline


def test_skeleton_columns_are_labelled_descriptive_exact_edge(tmp_path):
    result = run(tmp_path)
    assert "descriptive (skeleton exact-edge)" in result.stdout
    assert "undecidable" in result.stdout
