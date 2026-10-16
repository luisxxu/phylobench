import argparse
from pathlib import Path

from phylobench.regions import read_fasta
from phylobench.similarity import similarity, similarity_bin


def main():
    parser = argparse.ArgumentParser(
        description="Per-query similarity to the backbone from aligned queries."
    )
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--backbone", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    table = similarity(dict(read_fasta(args.queries)), dict(read_fasta(args.backbone)))
    table["bin"] = similarity_bin(table["matches"], table["n_columns"])
    table[["query", "sim", "bin"]].to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
