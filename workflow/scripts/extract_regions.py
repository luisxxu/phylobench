import argparse
from pathlib import Path

from phylobench.regions import (
    COLUMN_REGIONS,
    check_ungaps_to,
    column_window,
    cut_column_window,
    exclusion_counts,
    read_fasta,
    split_spans,
)


def main():
    parser = argparse.ArgumentParser(
        description="Cut one amplicon region from full-length queries."
    )
    parser.add_argument("--region", required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument(
        "--aligned",
        type=Path,
        help="shipped aligned queries; column-window regions only",
    )
    parser.add_argument("--out-fasta", type=Path, required=True)
    parser.add_argument("--out-exclusions", type=Path, required=True)
    parser.add_argument("--out-counts", type=Path, required=True)
    args = parser.parse_args()

    pairs = list(read_fasta(args.queries))
    if args.region in COLUMN_REGIONS:
        aligned = dict(read_fasta(args.aligned))
        check_ungaps_to(aligned, dict(pairs))
        window = column_window(aligned, args.region)
        print(f"{args.region} window: {window}")
        kept, exclusions = cut_column_window(aligned, window)
    else:
        spans, exclusions = split_spans(pairs, args.region)
        kept = {
            name: seq.upper()[slice(*spans[name])]
            for name, seq in pairs
            if name in spans
        }
    with args.out_fasta.open("w") as out:
        for name, _ in pairs:
            if name in kept:
                out.write(f">{name}\n{kept[name]}\n")
    args.out_exclusions.write_text(
        "".join(f"{name}\t{reason}\n" for name, reason in exclusions.items())
    )
    counts = {"extracted": len(kept), **exclusion_counts(exclusions)}
    args.out_counts.write_text(
        "reason\tcount\n" + "".join(f"{reason}\t{n}\n" for reason, n in counts.items())
    )


if __name__ == "__main__":
    main()
