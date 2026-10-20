import argparse
from pathlib import Path

import h5py
import treeswift

from phylobench.emp import asv_failures, asv_names, query_set

BIOM_FORMAT_VERSION = [2, 1]


def read_observation_ids(biom: Path) -> list[str]:
    with h5py.File(biom, "r") as f:
        version = f.attrs["format-version"].tolist()
        if version != BIOM_FORMAT_VERSION:
            raise SystemExit(f"{biom}: format-version {version}, expected 2.1")
        return [i.decode("ascii") for i in f["observation/ids"][:]]


def main():
    parser = argparse.ArgumentParser(
        description="Write a seeded, named subset of EMP ASVs as FASTA."
    )
    parser.add_argument("--biom", type=Path, required=True)
    parser.add_argument("--backbone-tree", type=Path, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    ids = read_observation_ids(args.biom)
    bad = asv_failures(ids)
    if bad:
        raise SystemExit(
            f"{len(bad)} observation IDs are not unique 150-bp ACGT: {bad[:3]}"
        )
    names = asv_names(args.n)
    leaves = set(
        treeswift.read_tree_newick(str(args.backbone_tree)).labels(internal=False)
    )
    collisions = sorted(leaves.intersection(names))
    if collisions:
        raise SystemExit(f"ASV names collide with backbone leaves: {collisions[:3]}")

    asvs = query_set(ids, args.n, args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        "".join(f">{name}\n{asv}\n" for name, asv in zip(names, asvs, strict=True))
    )


if __name__ == "__main__":
    main()
