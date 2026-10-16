import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

from phylobench.regions import read_fasta
from phylobench.tables import (
    RAW,
    load_arms,
    plan_cells,
    smoke_failures,
    width_failures,
)

ARMS_DIR = Path(__file__).parents[2] / "config" / "arms"


def main():
    parser = argparse.ArgumentParser(description="Assert the smoke run's acceptance.")
    parser.add_argument("--configfile", type=Path, required=True)
    parser.add_argument("--arm", required=True)
    args = parser.parse_args()

    config = yaml.safe_load(args.configfile.read_text())
    arms = load_arms(ARMS_DIR, only=args.arm)
    dataset, backbone = config["datasets"][0], config["backbones"][0]
    width = len(
        next(read_fasta(f"data/{dataset}/{backbone}/aligned-dna-sequences.fasta"))[1]
    )
    alns = {"hmmalign"} | {c.aln for c in plan_cells(config, arms) if c.aln != RAW}

    failures = []
    for region in config["regions"]:
        for aln in sorted(alns):
            path = Path(f"results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta")
            records = list(read_fasta(path)) if path.exists() else None
            failures += width_failures(str(path), records, width)

    queries = [n for n, _ in read_fasta(f"data/{dataset}/queries.fasta")]
    runs = pd.read_csv("results/runs.tsv", sep="\t")
    accuracy = pd.read_csv("results/accuracy.tsv", sep="\t")
    failures += smoke_failures(runs, accuracy, config, arms, queries)

    print("\n".join(failures) or "smoke ok")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
