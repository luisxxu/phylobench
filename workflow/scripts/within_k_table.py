import argparse
from pathlib import Path

import pandas as pd
import yaml

from phylobench.summaries import WITHIN_K, within_k_difference
from phylobench.tables import load_arms, plan_cells

ROOT = Path(__file__).parents[2]


def main():
    parser = argparse.ArgumentParser(
        description="Print the within-k non-inferiority table of a gg2-wol stage."
    )
    parser.add_argument("configs", type=Path, nargs="+")
    parser.add_argument("--accuracy", type=Path, default="results/accuracy.tsv")
    parser.add_argument("--heldout", type=Path, default="data/gg2-wol/heldout.tsv")
    parser.add_argument("--dataset", default="gg2-wol")
    parser.add_argument("-k", type=int, default=WITHIN_K)
    args = parser.parse_args()

    config = {}
    for path in args.configs:
        config.update(yaml.safe_load(path.read_text()))
    planned = [
        tuple(c)
        for c in plan_cells(config, load_arms(ROOT / "config" / "arms"))
        if c.dataset == args.dataset
    ]

    accuracy = pd.read_csv(
        args.accuracy, sep="\t", usecols=["run_id", "query", "edge_dist"]
    )
    accuracy = accuracy[accuracy["run_id"].str.startswith(f"{args.dataset}/")]
    heldout = pd.read_csv(args.heldout, sep="\t")
    groups = heldout[heldout["query"] == 1].set_index("leaf")["group"]

    table = within_k_difference(accuracy, groups, k=args.k, planned_cells=planned)
    print(table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
