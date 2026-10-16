import argparse
from pathlib import Path

from phylobench.regions import check_ungaps_to, project, read_fasta, split_spans


def main():
    parser = argparse.ArgumentParser(
        description="Oracle: the region's own GG2 alignment columns, gaps elsewhere."
    )
    parser.add_argument("--region", required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--aligned", type=Path, required=True)
    parser.add_argument("--backbone", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    full_length = dict(read_fasta(args.queries))
    aligned = dict(read_fasta(args.aligned))
    check_ungaps_to(aligned, full_length)
    width = len(next(read_fasta(args.backbone))[1])
    spans, _ = split_spans(full_length.items(), args.region)
    with args.out.open("w") as out:
        for name in full_length:
            if name not in spans:
                continue
            row = aligned[name]
            if len(row) != width:
                raise ValueError(
                    f"{name}: aligned width {len(row)} != backbone width {width}"
                )
            first, last = project(spans[name], row)
            out.write(
                f">{name}\n{'-' * first}{row[first:last]}{'-' * (width - last)}\n"
            )


if __name__ == "__main__":
    main()
