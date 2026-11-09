import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "within_k_table.py"
QUERIES = [f"G{i:09d}" for i in range(8)]


def accuracy_rows(arm, aln, edge_dist, dataset="gg2-wol"):
    return pd.DataFrame(
        {
            "run_id": f"{dataset}/B10k/full-length/{aln}/{arm}",
            "query": QUERIES,
            "edge_dist": edge_dist,
        }
    )


def test_prints_one_row_per_planned_challenger_cell(tmp_path):
    accuracy = pd.concat(
        [
            accuracy_rows("sepp-prod", "raw", [0.0] * 8),
            accuracy_rows("bscampp", "hmmalign", [0.0] * 6 + [9.0, 9.0]),
            accuracy_rows("sepp-prod", "raw", [9.0] * 8, dataset="gg2"),
        ]
    )
    accuracy.to_csv(tmp_path / "accuracy.tsv", sep="\t", index=False)
    heldout = pd.DataFrame(
        {
            "leaf": QUERIES,
            "group": ["clade1"] * 4 + [f"singleton{i}" for i in range(4)],
            "query": 1,
        }
    )
    heldout.to_csv(tmp_path / "heldout.tsv", sep="\t", index=False)
    config = tmp_path / "t28.yaml"
    config.write_text(
        "datasets: [gg2-wol]\nbackbones: [B10k]\nregions: [full-length]\n"
        "alns: [hmmalign]\nexclude:\n  - {arm: [pplacer, epa-ng, apples2,"
        " dipper-r, dipper-m, dipper-m-gpu, pplacer-scampp, epa-ng-scampp,"
        " sepp-pinned]}\n"
    )
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), str(config),
            "--accuracy", str(tmp_path / "accuracy.tsv"),
            "--heldout", str(tmp_path / "heldout.tsv"),
        ],
        capture_output=True,
        text=True,
    )  # fmt: skip
    assert result.returncode == 0, result.stderr
    (header, row) = result.stdout.strip().splitlines()
    cols = dict(zip(header.split(), row.split()))
    assert (cols["arm"], cols["n"], cols["n_clusters"]) == ("bscampp", "8", "5")
    assert cols["diff"] == "-0.2500"
    assert cols["verdict"] == "not-shown"
