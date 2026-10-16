import argparse
from pathlib import Path

import numpy as np

from phylobench.heldout import prune_alignment, prune_tree
from phylobench.simulate import (
    N_SITES,
    SEED,
    choose_heldout,
    embed_primer_sites,
    format_fasta,
    format_newick,
    format_raxml_info,
    format_truth,
    random_tree,
    simulate_alignment,
)

N_LEAVES = 10_000
CLADE_SIZES = (5, 10, 15, 20, 25) * 7
N_SINGLETONS = 500


def main():
    parser = argparse.ArgumentParser(
        description="Simulate the 10k-leaf true-tree panel (uDance-free cross-check)."
    )
    parser.add_argument("out_dir", type=Path, nargs="?", default=Path("data/sim10k"))
    out_dir = parser.parse_args().out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    reference = random_tree(N_LEAVES, rng)
    seqs = embed_primer_sites(simulate_alignment(reference, N_SITES, rng))
    groups = choose_heldout(reference, CLADE_SIZES, N_SINGLETONS, rng)
    query_names = [m for _, _, members in groups for m in members]
    held = set(query_names)
    backbone_seqs = prune_alignment(seqs, held)
    query_seqs = {name: seqs[name] for name in query_names}
    n_patterns = len(set(zip(*backbone_seqs.values(), strict=True)))

    files = {
        "reference_tree.nwk": format_newick(reference),
        "tree.nwk": format_newick(prune_tree(reference, held)),
        "aligned-dna-sequences.fasta": format_fasta(backbone_seqs),
        "queries.fasta": format_fasta(query_seqs),
        "queries.aligned.fasta": format_fasta(query_seqs),
        "raxml-info.txt": format_raxml_info(n_patterns),
        "truth.tsv": format_truth(groups),
    }
    for name, text in files.items():
        (out_dir / name).write_text(text)


if __name__ == "__main__":
    main()
