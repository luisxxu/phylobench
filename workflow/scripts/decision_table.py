import argparse
from pathlib import Path

import pandas as pd
import yaml

from phylobench.summaries import (
    WITHIN_K,
    decision_table,
    headline,
    verdict,
    within_k_difference,
)
from phylobench.tables import load_arms, plan_cells

ROOT = Path(__file__).parents[2]
CORE = {("raw", "sepp-pinned"), ("hmmalign", "bscampp")}


def main():
    parser = argparse.ArgumentParser(
        description="Print the T30 decision table and per-arm headline."
    )
    parser.add_argument("configs", type=Path, nargs="+")
    parser.add_argument("--accuracy", type=Path, default="results/accuracy.tsv")
    parser.add_argument("--skeleton", type=Path, default="results/skeleton.tsv")
    parser.add_argument("--heldout", type=Path, default="data/gg2-wol/heldout.tsv")
    parser.add_argument("--dataset", default="gg2-wol")
    args = parser.parse_args()

    config = {}
    for path in args.configs:
        config.update(yaml.safe_load(path.read_text()))
    planned = [
        tuple(c)
        for c in plan_cells(config, load_arms(ROOT / "config" / "arms"))
        if c.dataset == args.dataset
    ]

    def read(path):
        table = pd.read_csv(path, sep="\t", usecols=["run_id", "query", "edge_dist"])
        return table[table["run_id"].str.startswith(f"{args.dataset}/")]

    accuracy = read(args.accuracy)
    heldout = pd.read_csv(args.heldout, sep="\t")
    groups = heldout[heldout["query"] == 1].set_index("leaf")["group"]

    within = within_k_difference(accuracy, groups, k=WITHIN_K, planned_cells=planned)
    skeleton = within_k_difference(
        read(args.skeleton), groups, k=0, planned_cells=planned
    )
    decision = decision_table(within, verdict(accuracy), skeleton, CORE)

    print(
        "gate_025: the 0.25-edge mean gate is undecidable at these sample sizes.\n"
        "skel_diff, skel_lower: descriptive (skeleton exact-edge), k=0, no verdict.\n"
        "COST rows carry no verdict."
    )
    print(decision.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\nheadline")
    print(headline(decision).to_string(index=False))


if __name__ == "__main__":
    main()
