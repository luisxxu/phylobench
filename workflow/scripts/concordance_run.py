import argparse
import json
from pathlib import Path

import pandas as pd

from phylobench.concordance import score_concordance
from phylobench.regions import read_fasta
from phylobench.tables import CONCORDANCE_COLUMNS, concordance_rows


def is_ok(run_json: Path) -> bool:
    return json.loads(run_json.read_text())["status"] == "ok"


def read_jplace(jplace: Path) -> dict:
    text = jplace.read_text()
    if not text.strip():
        raise SystemExit(f"{jplace}: empty jplace from an ok placement run")
    return json.loads(text)


def main():
    parser = argparse.ArgumentParser(
        description="Edge distance between one cell's placements and sepp-prod's."
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--reference-run-id", required=True)
    parser.add_argument("--gate-run-json", type=Path, required=True)
    parser.add_argument("--reference-gate-run-json", type=Path, required=True)
    parser.add_argument("--jplace", type=Path, required=True)
    parser.add_argument("--reference-jplace", type=Path, required=True)
    parser.add_argument("--backbone-tree", type=Path, required=True)
    parser.add_argument("--region-queries", type=Path, required=True)
    parser.add_argument("--similarity", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    # A DNF'd placement phase left an empty jplace; its run.json carries the failure.
    if is_ok(args.gate_run_json) and is_ok(args.reference_gate_run_json):
        scores = score_concordance(
            read_jplace(args.jplace),
            read_jplace(args.reference_jplace),
            args.backbone_tree.read_text(),
            [name for name, _ in read_fasta(args.region_queries)],
        )
        sims = pd.read_csv(args.similarity, sep="\t")
        rows = concordance_rows(args.run_id, args.reference_run_id, scores, sims)
    pd.DataFrame(rows, columns=list(CONCORDANCE_COLUMNS)).to_csv(
        args.out, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
