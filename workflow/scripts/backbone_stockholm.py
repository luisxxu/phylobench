import argparse
from pathlib import Path

from phylobench.regions import read_fasta, stockholm_lines


def main():
    parser = argparse.ArgumentParser(
        description="Convert a backbone alignment to Stockholm with an all-x RF line."
    )
    parser.add_argument("--aligned", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    with args.out.open("w") as out:
        for line in stockholm_lines(read_fasta(args.aligned)):
            out.write(line + "\n")


if __name__ == "__main__":
    main()
