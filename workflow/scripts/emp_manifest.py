import argparse
from pathlib import Path

import h5py
import pandas as pd

from phylobench.emp import EMP_REGION
from phylobench.gg2 import MANIFEST_COLUMNS, bin_rows, file_row
from phylobench.measure import sha256_file
from phylobench.regions import read_fasta

BIOM = "data/emp/emp_deblur_150bp.release1.biom"
GG2_MANIFEST = "data/gg2/MANIFEST.tsv"
BACKBONE_FILES = ("tree.nwk", "aligned-dna-sequences.fasta", "raxml-info.txt")


def emp_file_row(path: str, n, backbone: str | None = None, **fields) -> dict:
    return {**file_row(path, backbone, n, sha256_file(Path(path))), **fields}


def backbone_rows(dataset: str, backbone: str, gg2: pd.DataFrame) -> list[dict]:
    known = gg2[gg2["kind"] == "file"].set_index("path")
    rows = []
    for name in BACKBONE_FILES:
        path = f"data/{dataset}/{backbone}/{name}"
        source = f"data/gg2/{backbone}/{name}"
        if source not in known.index:
            raise SystemExit(f"{GG2_MANIFEST} has no file row for {source}")
        n = known.at[source, "n"] or None
        row = emp_file_row(path, n, backbone, dataset=dataset)
        if known.at[source, "sha256"] != row["sha256"]:
            raise SystemExit(f"{path}: sha256 differs from {source} in {GG2_MANIFEST}")
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description="Write data/emp/MANIFEST.tsv.")
    parser.add_argument("--datasets", nargs="+", required=True)
    parser.add_argument("--backbones", nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    with h5py.File(BIOM, "r") as f:
        ids = f["observation/ids"].shape[0]
    gg2 = pd.read_csv(GG2_MANIFEST, sep="\t", dtype=str, keep_default_na=False)
    rows = [emp_file_row(BIOM, ids)]
    for dataset in args.datasets:
        queries = f"results/{dataset}/regions/{EMP_REGION}/queries.fasta"
        rows.append(
            emp_file_row(
                queries,
                sum(1 for _ in read_fasta(queries)),
                dataset=dataset,
                region=EMP_REGION,
            )
        )
        for backbone in args.backbones:
            rows.extend(backbone_rows(dataset, backbone, gg2))
            similarity = pd.read_csv(
                f"results/{dataset}/{backbone}/{EMP_REGION}/similarity.tsv", sep="\t"
            )
            rows += [
                {"dataset": dataset, **row}
                for row in bin_rows(backbone, EMP_REGION, similarity["bin"])
            ]
    manifest = pd.DataFrame(rows, columns=["dataset", *MANIFEST_COLUMNS])
    manifest.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
