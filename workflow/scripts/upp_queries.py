import argparse
import json
from pathlib import Path

from phylobench.regions import read_fasta, upp_query_rows


def main():
    parser = argparse.ArgumentParser(
        description="Validate UPP's masked alignment and keep the query rows."
    )
    parser.add_argument("--masked", type=Path, required=True)
    parser.add_argument("--backbone", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    # A DNF'd UPP left an empty masked file; the run.json row carries the failure.
    if json.loads(args.run_json.read_text())["status"] != "ok":
        args.out.write_text("")
        return
    queries = [name for name, _ in read_fasta(args.queries)]
    rows = upp_query_rows(
        dict(read_fasta(args.masked)), dict(read_fasta(args.backbone)), queries
    )
    print(f"{len(rows)} of {len(queries)} queries aligned by UPP")
    with args.out.open("w") as out:
        for name, row in rows:
            out.write(f">{name}\n{row}\n")


if __name__ == "__main__":
    main()
