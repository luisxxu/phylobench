import argparse
import json
from pathlib import Path

import pandas as pd

from phylobench.regions import read_fasta
from phylobench.tables import ACCURACY_COLUMNS, RUNS_COLUMNS, RunSpec, run_row


def query_count(spec: RunSpec) -> int | None:
    if not spec.region:
        return None
    path = f"results/{spec.dataset}/regions/{spec.region}/queries.fasta"
    return sum(1 for _ in read_fasta(path))


def main():
    parser = argparse.ArgumentParser(description="Collect run.json files into TSVs.")
    parser.add_argument("--specs", required=True, help="JSON list of RunSpec records")
    parser.add_argument("--out-runs", type=Path, required=True)
    parser.add_argument("--out-accuracy", type=Path, required=True)
    args = parser.parse_args()

    specs = [RunSpec(**record) for record in json.loads(args.specs)]
    runs = [
        run_row(spec, json.loads(Path(spec.run_json).read_text()), query_count(spec))
        for spec in specs
    ]
    pd.DataFrame(runs, columns=list(RUNS_COLUMNS)).to_csv(
        args.out_runs, sep="\t", index=False
    )
    scored = [
        pd.read_csv(Path(spec.run_json).with_name("accuracy.tsv"), sep="\t")
        for spec in specs
        if spec.phase == "evaluate"
    ]
    scored = [table for table in scored if len(table)]
    accuracy = (
        pd.concat(scored, ignore_index=True)
        if scored
        else pd.DataFrame(columns=list(ACCURACY_COLUMNS))
    )
    accuracy.to_csv(args.out_accuracy, sep="\t", index=False)


if __name__ == "__main__":
    main()
