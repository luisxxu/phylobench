import argparse
from pathlib import Path

import treeswift

from phylobench.gg2 import WOL_N_QUERIES, parse_lineages, wol_n_queries
from phylobench.regions import PRIMERS, read_fasta
from phylobench.wol import (
    clade_capacity,
    dry_run_capacity,
    genome_leaves,
    label_counts,
    length_quantiles,
    lineage_depths,
    neighbour_leaves,
    region_report,
    verdicts,
)


def leaf_labels(path: Path) -> list[str]:
    return [
        leaf.label
        for leaf in treeswift.read_tree_newick(path.read_text()).traverse_leaves()
    ]


def main():
    parser = argparse.ArgumentParser(description="Read-only S0 probe of the WoL set.")
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path, required=True)
    parser.add_argument("--backbone-dirs", type=Path, nargs="*", default=[])
    parser.add_argument("--n-queries", type=int, default=None)
    parser.add_argument("--neighbour-sample", type=int, default=2000)
    args = parser.parse_args()

    tree = treeswift.read_tree_newick((args.reference_dir / "tree.nwk").read_text())
    labels = [leaf.label for leaf in tree.traverse_leaves()]
    genomes = genome_leaves(labels)
    lineages = parse_lineages(args.taxonomy.read_text())
    counts = label_counts(labels)
    try:
        n_queries = args.n_queries or wol_n_queries(len(genomes))
    except ValueError:
        n_queries = WOL_N_QUERIES  # G1 prints STOP below
    print("P1-P3", counts)
    for directory in args.backbone_dirs:
        kept = genome_leaves(leaf_labels(directory / "tree.nwk"))
        print("P4", directory.name, "G leaves", len(kept))
    print("P5 G lineage depth", lineage_depths(genomes, lineages))
    for rank, (n, capacity) in clade_capacity(tree, lineages, genomes).items():
        print(f"P6 rank {rank}: {n} clades, query capacity {capacity}")

    seqs = dict(read_fasta(args.reference_dir / "aligned-dna-sequences.fasta"))
    rows = {name: seqs[name] for name in genomes}
    print("P7 G ungapped length", length_quantiles(rows))
    regions = {region: region_report(rows, region) for region in PRIMERS}
    for region, report in regions.items():
        print("P8", region, report)
    capacity = dry_run_capacity(tree, lineages, genomes, n_queries)
    print("P9", capacity)
    hit, sampled = neighbour_leaves(seqs, genomes, args.neighbour_sample)
    print(f"P10 {hit} of {sampled} sampled G leaves have a >=99.9 neighbour")
    for rule, verdict in verdicts(counts, capacity, regions, n_queries).items():
        print(rule, verdict)


if __name__ == "__main__":
    main()
