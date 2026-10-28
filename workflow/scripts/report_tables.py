import argparse
import os
import sys
from pathlib import Path

import pandas as pd

from phylobench.report import build_blocks, paste_blocks


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results", type=Path, default=os.environ.get("PHYLOBENCH_RESULTS", "results")
    )
    parser.add_argument("--licenses", type=Path, default="LICENSES.md")
    parser.add_argument(
        "--paste",
        type=Path,
        help="rewrite the table blocks of this results.md in place",
    )
    args = parser.parse_args()

    concordance_path = args.results / "concordance.tsv"
    concordance = (
        read_tsv(concordance_path, usecols=["run_id", "edge_dist"])
        if concordance_path.exists()
        else None
    )
    blocks = build_blocks(
        read_tsv(args.results / "runs.tsv"),
        read_tsv(args.results / "accuracy.tsv"),
        concordance,
        args.licenses.read_text(),
    )
    if args.paste:
        args.paste.write_text(paste_blocks(args.paste.read_text(), blocks))
    else:
        sys.stdout.write("\n\n".join(blocks.values()) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
