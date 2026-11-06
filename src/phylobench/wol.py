"""S0 feasibility probe of the WoL genome held-out set (docs/runbook-t31.md).

Pure and read-only: each function turns the shipped tree, alignment and taxonomy into
numbers the pre-registered stop/go rules G1-G6 are read from.
"""

import re
from collections import Counter
from collections.abc import Mapping
from statistics import median

import numpy as np
import treeswift

from phylobench.gg2 import (
    MAX_LEAVES,
    MIN_LEAVES,
    PER_CLADE,
    WOL_LEAF,
    WOL_MIN_ELIGIBLE,
    WOL_RANKS,
    WOL_SEED,
    near_identical,
    wol_n_queries,
)
from phylobench.heldout import choose_groups, monophyletic_taxa, sample_queries
from phylobench.regions import (
    COLUMN_SUFFIX,
    amplicon_span,
    column_window,
    cut_column_window,
    ungap,
)

_G_ANYWHERE = re.compile(r"G\d{9}")
_ACCESSION = re.compile(r"(?:RS|GB)_GC[AF]_(\d{9})")
MIN_C_ROWS = 30


def genome_leaves(labels) -> frozenset[str]:
    return frozenset(label for label in labels if WOL_LEAF.match(label))


def label_counts(labels: list[str]) -> dict:
    """P1-P3: totals, G leaves, G-like labels that fail the pattern, RS_/GB_ twins."""
    genomes = genome_leaves(labels)
    odd = Counter(
        _G_ANYWHERE.sub("G#", label)
        for label in labels
        if _G_ANYWHERE.search(label) and label not in genomes
    )
    accessions = {m[1] for label in labels if (m := _ACCESSION.match(label))}
    return {
        "leaves": len(labels),
        "G": len(genomes),
        "G_like_not_matching": odd.most_common(5),
        "G_with_RS_GB_twin": sum(label[1:] in accessions for label in genomes),
    }


def lineage_depths(genomes, lineages: Mapping[str, tuple[str, ...]]) -> list:
    """P5: (lineage depth, G leaves) pairs."""
    return sorted(Counter(len(lineages.get(g, ())) for g in genomes).items())


def clade_capacity(
    tree: treeswift.Tree,
    lineages: Mapping[str, tuple[str, ...]],
    genomes: frozenset[str],
    ranks=WOL_RANKS,
) -> dict[int, tuple[int, int]]:
    """P6: per rank, (taxa of 25-2,000 leaves with >=5 G, their query capacity)."""
    capacity = {}
    for rank in sorted(ranks):
        counts = []
        for node in monophyletic_taxa(tree, dict(lineages), rank).values():
            members = [leaf.label for leaf in node.traverse_leaves()]
            n_genomes = sum(m in genomes for m in members)
            if MIN_LEAVES <= len(members) <= MAX_LEAVES and (
                n_genomes >= WOL_MIN_ELIGIBLE
            ):
                counts.append(min(n_genomes, PER_CLADE))
        capacity[rank] = (len(counts), sum(counts))
    return capacity


def length_quantiles(rows: Mapping[str, str]) -> dict[str, float]:
    """P7: ungapped lengths."""
    lengths = np.array([len(ungap(row)) for row in rows.values()])
    q10, q90 = np.quantile(lengths, [0.1, 0.9])
    return {
        "min": lengths.min(), "q10": q10, "median": median(lengths.tolist()),
        "q90": q90, "max": lengths.max(),
    }  # fmt: skip


def region_report(rows: Mapping[str, str], region: str) -> dict:
    """P8: primer-hit counts, the `c` window and rows it keeps, and the G5 verdict."""
    hits = Counter(
        "ok" if not isinstance(span := amplicon_span(ungap(row), region), str) else span
        for row in rows.values()
    )
    try:
        window = column_window(rows, region + COLUMN_SUFFIX)
        kept = len(cut_column_window(rows, window)[0])
    except ValueError:
        window, kept = None, 0
    return {
        "primer_hits": dict(hits),
        "window": window,
        "c_rows_kept": kept,
        "c_dropped_g5": hits["ok"] < MIN_C_ROWS,
    }


def dry_run_capacity(
    tree: treeswift.Tree,
    lineages: Mapping[str, tuple[str, ...]],
    genomes: frozenset[str],
    n_queries: int,
    ranks=WOL_RANKS,
) -> dict:
    """P9: draw the clades and singletons of the build; report what could be drawn.

    `g_outside_clades` bounds the singletons; `singletons_ok` is whether `n_queries`
    minus the clade queries could really be drawn disjointly (G3).
    """
    rng = np.random.default_rng(WOL_SEED)
    common = {
        "max_leaves": MAX_LEAVES, "min_leaves": MIN_LEAVES,
        "eligible": genomes, "min_eligible": WOL_MIN_ELIGIBLE,
    }  # fmt: skip
    try:
        clades = choose_groups(tree, dict(lineages), ranks, 0, rng=rng, **common)
    except ValueError as error:
        return {"clades_ok": False, "error": str(error)}
    clade_queries = len(sample_queries(clades, rng, PER_CLADE, genomes))
    in_clades = {m for g in clades for m in g.members}
    wanted = max(n_queries - clade_queries, 0)
    rng = np.random.default_rng(WOL_SEED)
    try:
        choose_groups(tree, dict(lineages), ranks, wanted, rng=rng, **common)
        singletons_ok, error = True, ""
    except ValueError as exc:
        singletons_ok, error = False, str(exc)
    return {
        "clades_ok": True,
        "clades": len(clades),
        "clade_queries": clade_queries,
        "singletons_wanted": wanted,
        "g_outside_clades": len(genomes - in_clades),
        "singletons_ok": singletons_ok,
        "error": error,
    }


def neighbour_leaves(
    seqs: Mapping[str, str], genomes: frozenset[str], sample: int, seed: int = 0
) -> tuple[int, int]:
    """P10: of `sample` random G leaves, how many have a >=99.9 neighbour among all
    other leaves (the sample itself excluded, as when seeding)."""
    names = sorted(genomes & seqs.keys())
    picks = np.random.default_rng(seed).permutation(len(names))[:sample]
    chosen = {names[i] for i in picks}
    others = {n: s for n, s in seqs.items() if n not in chosen}
    found = near_identical({n: seqs[n] for n in chosen}, others)
    return int(sum(found)), len(chosen)


MIN_V4_QUERIES = 2000
MIN_CLADE_QUERIES = 500


def verdicts(
    counts: dict, capacity: dict, regions: dict[str, dict], n_queries: int
) -> dict[str, str]:
    """Pre-registered G1-G5 (docs/runbook-t31.md) from the probe numbers.

    G2 cannot be settled by a count alone, so copy-suffixed labels give REVIEW.
    """
    try:
        g1 = f"GO N={wol_n_queries(counts['G'])}"
    except ValueError as error:
        g1 = str(error)
    g2 = (
        "GO"
        if not counts["G_like_not_matching"]
        else f"REVIEW {counts['G_like_not_matching']}"
    )
    if not capacity["clades_ok"]:
        g3 = f"STOP clades: {capacity['error']}"
    else:
        short = capacity["clade_queries"] < MIN_CLADE_QUERIES or (
            not capacity["singletons_ok"]
        )
        g3 = f"{'STOP' if short else 'GO'} {capacity}"
    v4 = regions["V4-150"]["primer_hits"].get("ok", 0)
    expected = v4 / counts["G"] * n_queries
    g4 = f"{'STOP' if expected < MIN_V4_QUERIES else 'GO'} expect {expected:.0f} V4-150"
    dropped = sorted(r for r, rep in regions.items() if rep["c_dropped_g5"])
    return {"G1": g1, "G2": g2, "G3": g3, "G4": g4, "G5": f"drop c twin: {dropped}"}
