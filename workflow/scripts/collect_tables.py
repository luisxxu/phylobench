import argparse
import json
from pathlib import Path

import pandas as pd

from phylobench.regions import read_fasta
from phylobench.tables import (
    ACCURACY_COLUMNS,
    CONCORDANCE_COLUMNS,
    CONCORDANCE_RUN,
    RUNS_COLUMNS,
    RunSpec,
    run_row,
    score_table,
)


def query_count(spec: RunSpec) -> int | None:
    if not spec.region:
        return None
    path = f"results/{spec.dataset}/regions/{spec.region}/queries.fasta"
    return sum(1 for _ in read_fasta(path))


def write_scores(specs: list[RunSpec], columns, out: Path) -> None:
    tables = [pd.read_csv(score_table(spec), sep="\t") for spec in specs]
    tables = [table for table in tables if len(table)]
    scores = (
        pd.concat(tables, ignore_index=True)
        if tables
        else pd.DataFrame(columns=list(columns))
    )
    scores.to_csv(out, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Collect run.json files into TSVs.")
    parser.add_argument(
        "--specs", type=Path, required=True, help="JSON file of RunSpec records"
    )
    parser.add_argument("--out-runs", type=Path, required=True)
    parser.add_argument("--out-accuracy", type=Path, required=True)
    parser.add_argument("--out-concordance", type=Path, required=True)
    args = parser.parse_args()

    specs = [RunSpec(**record) for record in json.loads(args.specs.read_text())]
    runs = [
        run_row(spec, json.loads(Path(spec.run_json).read_text()), query_count(spec))
        for spec in specs
    ]
    pd.DataFrame(runs, columns=list(RUNS_COLUMNS)).to_csv(
        args.out_runs, sep="\t", index=False
    )
    evaluated = [spec for spec in specs if spec.phase == "evaluate"]
    write_scores(
        [s for s in evaluated if Path(s.run_json).name != CONCORDANCE_RUN],
        ACCURACY_COLUMNS,
        args.out_accuracy,
    )
    write_scores(
        [s for s in evaluated if Path(s.run_json).name == CONCORDANCE_RUN],
        CONCORDANCE_COLUMNS,
        args.out_concordance,
    )


if __name__ == "__main__":
    main()
