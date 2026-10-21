import argparse
import sys
from pathlib import Path

import pandas as pd

from phylobench.tables import timing_summary


def main():
    parser = argparse.ArgumentParser(
        description="Place-phase timing summary over one or more runs.tsv files."
    )
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    runs = pd.concat([pd.read_csv(p, sep="\t") for p in args.runs], ignore_index=True)
    summary = timing_summary(runs)
    summary.to_csv(args.out or sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
