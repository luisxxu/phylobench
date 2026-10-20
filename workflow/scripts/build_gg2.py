import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd
import treeswift

from phylobench.gg2 import (
    B10K,
    HELDOUT_COLUMNS,
    MANIFEST_COLUMNS,
    TARGET,
    bin_rows,
    build_backbones,
    check_backbone,
    choose_seeded_groups,
    file_row,
    heldout_rows,
    parse_lineages,
    reference_tree,
    region_bins,
    similarity_memory_mb,
)
from phylobench.measure import sha256_file
from phylobench.regions import REGIONS, read_fasta, ungap
from phylobench.similarity import TOP_BIN
from phylobench.simulate import format_fasta, format_newick


def log(message):
    print(message, file=sys.stderr, flush=True)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def main():
    parser = argparse.ArgumentParser(description="Build the GG2 backbones and queries.")
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path, required=True)
    parser.add_argument("--downloads", type=Path, nargs="+", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mem-mb", type=int, required=True)
    parser.add_argument("--b10k", type=int, default=B10K)
    args = parser.parse_args()

    seqs = dict(read_fasta(args.reference_dir / "aligned-dna-sequences.fasta"))
    tree = treeswift.read_tree_newick((args.reference_dir / "tree.nwk").read_text())
    lineages = parse_lineages(args.taxonomy.read_text())
    width = len(next(iter(seqs.values())))
    needed_mb = similarity_memory_mb(len(seqs), width)
    log(f"R={len(seqs)} L={width} similarity needs about {needed_mb:.0f} MB")
    if needed_mb > args.mem_mb:
        sys.exit(f"need {needed_mb:.0f} MB for similarity, only {args.mem_mb} MB given")

    groups, queries, seeds = choose_seeded_groups(tree, seqs, lineages)
    held_out = {m for g in groups for m in g.members}
    log(f"held out {len(held_out)} leaves, {len(queries)} queries, {len(seeds)} seeded")

    backbones = build_backbones(tree, seqs, held_out, args.b10k)
    log(f"removed {len(backbones.duplicates)} duplicate rows from Bfull")
    built = {
        "Bfull": (backbones.full_tree, backbones.full_seqs),
        "B10k": (backbones.small_tree, backbones.small_seqs),
    }
    for backbone_tree, backbone_seqs in built.values():
        check_backbone(backbone_tree, backbone_seqs, held_out)

    out = args.out_dir
    manifest = [
        file_row(str(path), None, None, sha256_file(path)) for path in args.downloads
    ]
    for name, (backbone_tree, backbone_seqs) in built.items():
        directory = out / name
        keep = {leaf.label for leaf in backbone_tree.traverse_leaves()} | set(queries)
        files = {
            "tree.nwk": format_newick(backbone_tree),
            "aligned-dna-sequences.fasta": format_fasta(backbone_seqs),
            "reference_tree.nwk": format_newick(reference_tree(tree, keep)),
        }
        for filename, text in files.items():
            write_text(directory / filename, text)
        shutil.copyfile(
            args.reference_dir / "raxml-info.txt", directory / "raxml-info.txt"
        )
        leaf_counts = {
            "tree.nwk": len(backbone_seqs),
            "aligned-dna-sequences.fasta": len(backbone_seqs),
            "reference_tree.nwk": len(keep),
            "raxml-info.txt": None,
        }
        for filename, n in leaf_counts.items():
            path = directory / filename
            manifest.append(file_row(str(path), name, n, sha256_file(path)))

    aligned = {q: seqs[q] for q in queries}
    ungapped = {q: ungap(row) for q, row in aligned.items()}
    write_text(out / "queries.fasta", format_fasta(ungapped))
    write_text(out / "queries.aligned.fasta", format_fasta(aligned))
    heldout = pd.DataFrame(
        heldout_rows(groups, queries, seeds), columns=HELDOUT_COLUMNS
    )
    heldout.to_csv(out / "heldout.tsv", sep="\t", index=False)
    for filename, n in (
        ("queries.fasta", len(queries)),
        ("queries.aligned.fasta", len(queries)),
        ("heldout.tsv", len(heldout)),
    ):
        manifest.append(
            file_row(str(out / filename), None, n, sha256_file(out / filename))
        )

    for name, (_, backbone_seqs) in built.items():
        for region in REGIONS:
            bins = region_bins(ungapped, aligned, backbone_seqs, region)
            manifest += bin_rows(name, region, bins)
            log(f"{name} {region}: {len(bins)} queries binned")
            top = bins.count(TOP_BIN)
            if region == "full-length" and top < TARGET:
                log(f"WARNING {name} full-length >=99.9 has {top} queries")

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(manifest, columns=MANIFEST_COLUMNS).to_csv(
        args.manifest, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
