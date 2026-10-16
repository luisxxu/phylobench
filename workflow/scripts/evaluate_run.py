import argparse
import json
from pathlib import Path

import pandas as pd

from phylobench.accuracy import score_jplace
from phylobench.regions import read_fasta
from phylobench.tables import ACCURACY_COLUMNS, accuracy_rows


def main():
    parser = argparse.ArgumentParser(
        description="Score one cell's jplace against the reference tree."
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--gate-run-json", type=Path, required=True)
    parser.add_argument("--jplace", type=Path, required=True)
    parser.add_argument("--backbone-tree", type=Path, required=True)
    parser.add_argument("--reference-tree", type=Path, required=True)
    parser.add_argument("--region-queries", type=Path, required=True)
    parser.add_argument("--similarity", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    # A DNF'd placement phase left an empty jplace; its run.json carries the failure.
    if json.loads(args.gate_run_json.read_text())["status"] == "ok":
        text = args.jplace.read_text()
        if not text.strip():
            raise SystemExit(f"{args.jplace}: empty jplace from an ok placement run")
        queries = [name for name, _ in read_fasta(args.region_queries)]
        scores = score_jplace(
            json.loads(text),
            args.backbone_tree.read_text(),
            args.reference_tree.read_text(),
            queries,
        )
        sims = pd.read_csv(args.similarity, sep="\t")
        rows = accuracy_rows(args.run_id, scores, sims, label="accuracy")
    pd.DataFrame(rows, columns=list(ACCURACY_COLUMNS)).to_csv(
        args.out, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
