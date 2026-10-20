"""GG2 2022.10 dataset construction (design.md section 4): pure functions.

Held-out clades and singletons come out of the shipped tree and alignment; Bfull* is
what remains (exact duplicate rows removed), B10k* a uniform subsample of Bfull*.
"""

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import treeswift

from phylobench.heldout import (
    HeldOutGroup,
    choose_groups,
    prune_alignment,
    prune_tree,
    sample_queries,
)
from phylobench.regions import project, split_spans
from phylobench.similarity import (
    BIN_LABELS,
    MIN_BIN_SIZE,
    TOP_BIN,
    similarity,
    similarity_bin,
)

RANKS = {1: 2, 2: 3, 3: 4, 4: 6, 5: 6}  # phylum, class, order, family, genus counts
MIN_LEAVES = 25
MAX_LEAVES = 2000
PER_CLADE = 25
N_QUERIES = 1000
SEED = 2022
N_CANDIDATES = 2000
TARGET = 100  # queries needed at >=99.9 before a bin stops being flagged
B10K = 10_000

MANIFEST_COLUMNS = (
    "kind", "backbone", "region", "bin", "path", "n", "sha256", "flagged"
)  # fmt: skip
HELDOUT_COLUMNS = ("leaf", "group", "group_type", "rank", "taxon", "query", "seeded")


def parse_lineages(tsv_text: str) -> dict[str, tuple[str, ...]]:
    """Feature ID -> ranks of a QIIME taxonomy.tsv, cut at the first empty rank.

    Taxonomy IDs spell the tree's underscores as dashes, so IDs are mapped back.
    """
    lineages = {}
    for line in tsv_text.splitlines()[1:]:
        feature, taxon = line.split("\t")
        ranks = []
        for rank in taxon.split("; "):
            if rank.endswith("__"):
                break
            ranks.append(rank)
        lineages[feature.replace("-", "_")] = tuple(ranks)
    return lineages


def _row_key(seq: str) -> str:
    return seq.upper().replace(".", "-")


def duplicate_rows(seqs: Mapping[str, str]) -> set[str]:
    """Names of rows equal to an earlier row in iteration order."""
    seen: set[str] = set()
    duplicates: set[str] = set()
    for name, seq in seqs.items():
        key = _row_key(seq)
        if key in seen:
            duplicates.add(name)
        seen.add(key)
    return duplicates


def check_backbone(
    tree: treeswift.Tree, seqs: Mapping[str, str], heldout: set[str]
) -> None:
    """Raise ValueError unless the backbone is what FastTree and every arm expect."""
    leaves = {leaf.label for leaf in tree.traverse_leaves()}
    if leaves != seqs.keys():
        raise ValueError("tree leaves differ from alignment rows")
    if not leaves.isdisjoint(heldout):
        raise ValueError("held-out leaves present in the backbone")
    if len({len(seq) for seq in seqs.values()}) > 1:
        raise ValueError("alignment rows have unequal widths")
    if duplicate_rows(seqs):
        raise ValueError("duplicate alignment rows")
    if len(tree.root.children) not in (2, 3):
        raise ValueError("root is not bifurcating or trifurcating")
    for node in tree.traverse_preorder(leaves=False):
        if node is not tree.root and len(node.children) != 2:
            raise ValueError("tree is not strictly bifurcating")


def seeds_needed(hits: Sequence[int], n_available: int, target: int = TARGET) -> int:
    """Fewest trailing singletons to replace by seeds so hits reach `target`."""
    for k in range(len(hits) + 1):
        if sum(hits[: len(hits) - k]) + k >= target:
            if k > n_available:
                raise ValueError(
                    f"{k} seeds needed but only {n_available} candidates qualify; "
                    "raise n_candidates"
                )
            return k
    raise ValueError(f"cannot reach {target} hits with {len(hits)} singletons")


def similarity_memory_mb(n_refs: int, n_columns: int) -> float:
    return 16 * n_refs * n_columns / 2**20


def near_identical(queries: dict[str, str], refs: dict[str, str]) -> list[bool]:
    """Per query: full-length similarity to `refs` in the top (>=99.9) bin."""
    table = similarity(queries, refs)
    return list(similarity_bin(table["matches"], table["n_columns"]) == TOP_BIN)


def region_bins(
    ungapped: dict[str, str],
    aligned: dict[str, str],
    refs: dict[str, str],
    region: str,
) -> list[str]:
    """Bin of each query that has `region`, scored inside its own window."""
    spans, _ = split_spans(ungapped.items(), region)
    rows = {name: aligned[name] for name in spans}
    windows = {name: project(spans[name], rows[name]) for name in spans}
    table = similarity(rows, refs, windows)
    return list(similarity_bin(table["matches"], table["n_columns"]))


def _draw(
    tree, lineages, rank_counts, n_singletons, seed, *, max_leaves, min_leaves,
    per_clade, preferred=(), n_preferred=0,
):  # fmt: skip
    rng = np.random.default_rng(seed)
    groups = choose_groups(
        tree,
        lineages,
        rank_counts,
        n_singletons,
        max_leaves=max_leaves,
        rng=rng,
        min_leaves=min_leaves,
        preferred=preferred,
        n_preferred=n_preferred,
    )
    return groups, sample_queries(groups, rng, per_clade)


def choose_seeded_groups(
    tree: treeswift.Tree,
    seqs: dict[str, str],
    lineages: dict[str, tuple[str, ...]],
    *,
    rank_counts: dict[int, int] = RANKS,
    n_queries: int = N_QUERIES,
    max_leaves: int = MAX_LEAVES,
    min_leaves: int = MIN_LEAVES,
    per_clade: int = PER_CLADE,
    seed: int = SEED,
    n_candidates: int = N_CANDIDATES,
    target: int = TARGET,
) -> tuple[list[HeldOutGroup], list[str], set[str]]:
    """Held-out groups, the queries, and the leaves that were seeded.

    Singletons are drawn at random. If fewer than `target` of them have a >=99.9
    neighbour left in the backbone, the trailing singletons are replaced by
    candidates that do (design.md section 4 step 2).
    """
    shape = {
        "max_leaves": max_leaves,
        "min_leaves": min_leaves,
        "per_clade": per_clade,
    }
    clades, _ = _draw(tree, lineages, rank_counts, 0, seed, **shape)
    n_singletons = n_queries - sum(min(len(g.members), per_clade) for g in clades)
    if n_singletons < 0:
        raise ValueError("clade queries alone exceed n_queries")

    groups, queries = _draw(tree, lineages, rank_counts, n_singletons, seed, **shape)
    leaves = [g.members[0] for g in groups if g.group_type == "singleton"]
    held_out = {m for g in groups for m in g.members}
    remaining = prune_alignment(seqs, held_out)
    hits = [int(h) for h in near_identical({s: seqs[s] for s in leaves}, remaining)]
    if sum(hits) >= target:
        return groups, queries, set()

    names = list(remaining)
    picks = np.random.default_rng([seed, 1]).permutation(len(names))[:n_candidates]
    candidates = [names[i] for i in picks]
    chosen = set(candidates)
    others = {n: s for n, s in remaining.items() if n not in chosen}
    found = near_identical({c: seqs[c] for c in candidates}, others)
    seeds = [c for c, hit in zip(candidates, found, strict=True) if hit]

    k = seeds_needed(hits, len(seeds), target)
    groups, queries = _draw(
        tree, lineages, rank_counts, n_singletons, seed,
        preferred=seeds, n_preferred=k, **shape,
    )  # fmt: skip
    leaves = [g.members[0] for g in groups if g.group_type == "singleton"]
    return groups, queries, set(leaves[:k]) & set(seeds)


@dataclass(frozen=True)
class Backbones:
    full_tree: treeswift.Tree
    full_seqs: dict[str, str]
    small_tree: treeswift.Tree
    small_seqs: dict[str, str]
    duplicates: set[str]


def build_backbones(
    tree: treeswift.Tree,
    seqs: dict[str, str],
    held_out: set[str],
    b10k_size: int = B10K,
    seed: int = SEED,
) -> Backbones:
    """Bfull* (held-out and then duplicate rows removed) and B10k* (uniform subset)."""
    pruned_tree = prune_tree(tree, held_out)
    pruned_seqs = prune_alignment(seqs, held_out)
    order = [leaf.label for leaf in pruned_tree.traverse_leaves()]
    duplicates = duplicate_rows({name: pruned_seqs[name] for name in order})
    full_tree = prune_tree(pruned_tree, duplicates) if duplicates else pruned_tree
    full_seqs = prune_alignment(pruned_seqs, duplicates)

    names = [name for name in order if name not in duplicates]
    picks = np.random.default_rng(seed).choice(len(names), b10k_size, replace=False)
    keep = {names[i] for i in picks}
    small_tree = prune_tree(full_tree, set(names) - keep)
    small_seqs = {n: s for n, s in full_seqs.items() if n in keep}
    return Backbones(full_tree, full_seqs, small_tree, small_seqs, duplicates)


def reference_tree(tree: treeswift.Tree, keep: set[str]) -> treeswift.Tree:
    """The shipped tree restricted to `keep` (a backbone plus its queries)."""
    labels = {leaf.label for leaf in tree.traverse_leaves()}
    return prune_tree(tree, labels - keep)


def heldout_rows(
    groups: Iterable[HeldOutGroup], queries: Iterable[str], seeds: set[str]
) -> list[dict]:
    query_set = set(queries)
    return [
        {
            "leaf": leaf,
            "group": group.group_id,
            "group_type": group.group_type,
            "rank": "" if group.rank is None else group.rank,
            "taxon": group.taxon,
            "query": int(leaf in query_set),
            "seeded": int(leaf in seeds),
        }
        for group in groups
        for leaf in group.members
    ]


def bin_rows(backbone: str, region: str, bins: Iterable[str]) -> list[dict]:
    """One manifest row per similarity bin, empty bins included."""
    counts = Counter(bins)
    return [
        {
            "kind": "bin",
            "backbone": backbone,
            "region": region,
            "bin": label,
            "path": "",
            "n": counts[label],
            "sha256": "",
            "flagged": int(counts[label] < MIN_BIN_SIZE),
        }
        for label in (*BIN_LABELS, TOP_BIN)
    ]


def file_row(path: str, backbone: str | None, n: int | None, sha256: str) -> dict:
    return {
        "kind": "file",
        "backbone": backbone or "",
        "region": "",
        "bin": "",
        "path": path,
        "n": "" if n is None else n,
        "sha256": sha256,
        "flagged": "",
    }
