import argparse
import json
from pathlib import Path

from phylobench.regions import mask_inserts, read_fasta


def main():
    parser = argparse.ArgumentParser(
        description="Mask hmmalign output to the backbone's match columns."
    )
    parser.add_argument("--afa", type=Path, required=True)
    parser.add_argument("--backbone", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    # A DNF'd hmmalign left an empty afa; the run.json row carries the failure.
    if json.loads(args.run_json.read_text())["status"] != "ok":
        args.out.write_text("")
        return
    width = len(next(read_fasta(args.backbone))[1])
    with args.out.open("w") as out:
        for name, row in read_fasta(args.afa):
            masked = mask_inserts(row)
            if len(masked) != width:
                raise ValueError(
                    f"{name}: masked width {len(masked)} != backbone width {width}"
                )
            out.write(f">{name}\n{masked}\n")


if __name__ == "__main__":
    main()
