import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml
from snakemake.utils import validate

from phylobench.gates import stage_report
from phylobench.tables import load_arms

ROOT = Path(__file__).parents[2]


def main():
    parser = argparse.ArgumentParser(
        description="Gate a stage's results/*.tsv against the planned matrix."
    )
    parser.add_argument(
        "configs", type=Path, nargs="+", help="configfiles in snakemake merge order"
    )
    args = parser.parse_args()

    config = {}
    for path in args.configs:
        config.update(yaml.safe_load(path.read_text()))
    config["arms"] = load_arms(ROOT / "config" / "arms")
    validate(config, str(ROOT / "config" / "schema.yaml"))

    runs = pd.read_csv("results/runs.tsv", sep="\t")
    accuracy = pd.read_csv("results/accuracy.tsv", sep="\t", usecols=["run_id"])
    concordance = pd.read_csv("results/concordance.tsv", sep="\t", usecols=["run_id"])
    lines, failed = stage_report(runs, accuracy, concordance, config, config["arms"])
    print("\n".join(lines))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
