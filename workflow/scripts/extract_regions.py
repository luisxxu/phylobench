import argparse
from pathlib import Path

from phylobench.regions import exclusion_counts, read_fasta, split_spans


def main():
    parser = argparse.ArgumentParser(
        description="Cut one amplicon region from full-length queries."
    )
    parser.add_argument("--region", required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--out-fasta", type=Path, required=True)
    parser.add_argument("--out-exclusions", type=Path, required=True)
    parser.add_argument("--out-counts", type=Path, required=True)
    args = parser.parse_args()

    pairs = list(read_fasta(args.queries))
    spans, exclusions = split_spans(pairs, args.region)
    with args.out_fasta.open("w") as out:
        for name, seq in pairs:
            if name in spans:
                start, end = spans[name]
                out.write(f">{name}\n{seq.upper()[start:end]}\n")
    args.out_exclusions.write_text(
        "".join(f"{name}\t{reason}\n" for name, reason in exclusions.items())
    )
    counts = {"extracted": len(spans), **exclusion_counts(exclusions)}
    args.out_counts.write_text(
        "reason\tcount\n" + "".join(f"{reason}\t{n}\n" for reason, n in counts.items())
    )


if __name__ == "__main__":
    main()
