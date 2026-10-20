import hashlib
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import treeswift
from tree_helpers import coalescent_tree

from phylobench.accuracy import prepare_backbone, true_positions
from phylobench.gg2 import (
    MANIFEST_COLUMNS,
    bin_rows,
    build_backbones,
    check_backbone,
    choose_seeded_groups,
    duplicate_rows,
    file_row,
    heldout_rows,
    near_identical,
    parse_lineages,
    reference_tree,
    region_bins,
    seeds_needed,
)
from phylobench.heldout import choose_groups, prune_alignment, sample_queries
from phylobench.regions import read_fasta, ungap
from phylobench.similarity import BIN_LABELS, TOP_BIN

ROOT = Path(__file__).parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
FIXTURE = ROOT / "tests" / "fixtures" / "tiny"


def tree_of(newick):
    return treeswift.read_tree_newick(newick)


def test_parse_lineages_truncates_at_first_empty_rank_and_maps_dashes():
    text = (
        "Feature ID\tTaxon\n"
        "AB-1\td__Bacteria; p__Firm; c__Clos; o__; f__Fam; g__; s__\n"
        "CD2\td__Archaea; p__Eury; c__; o__; f__; g__; s__\n"
    )
    assert parse_lineages(text) == {
        "AB_1": ("d__Bacteria", "p__Firm", "c__Clos"),
        "CD2": ("d__Archaea", "p__Eury"),
    }


def test_duplicate_rows_keeps_first_and_ignores_case_and_dot_gaps():
    seqs = {"a": "AC-GT", "b": "ac.gt", "c": "ACGGT", "d": "AC-GT"}
    assert duplicate_rows(seqs) == {"b", "d"}


def good_backbone():
    tree = tree_of("((a:1,b:1):1,(c:1,d:1):1,e:1);")
    seqs = {"a": "AAAA", "b": "AAAC", "c": "AACC", "d": "ACCC", "e": "CCCC"}
    return tree, seqs


def test_check_backbone_accepts_valid():
    tree, seqs = good_backbone()
    check_backbone(tree, seqs, {"z"})


@pytest.mark.parametrize(
    ("mutate", "match"),
    [
        (lambda t, s, h: s.pop("e"), "leaves"),
        (lambda t, s, h: h.add("a"), "held-out"),
        (lambda t, s, h: s.update(b="AAAA"), "duplicate"),
        (lambda t, s, h: s.update(b="AAA"), "width"),
    ],
)
def test_check_backbone_failures(mutate, match):
    tree, seqs = good_backbone()
    held = {"z"}
    mutate(tree, seqs, held)
    with pytest.raises(ValueError, match=match):
        check_backbone(tree, seqs, held)


def test_check_backbone_rejects_polytomy_and_wide_root():
    seqs = {
        c: s
        for c, s in zip(
            "abcdef", ["AAAA", "AAAC", "AACC", "ACCC", "CCCC", "CCCA"], strict=True
        )
    }
    polytomy = tree_of("((a:1,b:1,c:1):1,(d:1,e:1):1,f:1);")
    with pytest.raises(ValueError, match="bifurcating"):
        check_backbone(polytomy, seqs, set())
    wide_root = tree_of("((a:1,b:1):1,(c:1,d:1):1,e:1,f:1);")
    with pytest.raises(ValueError, match="bifurcating"):
        check_backbone(wide_root, seqs, set())


def test_check_backbone_accepts_two_child_root():
    seqs = {"a": "AAAA", "b": "AAAC", "c": "AACC"}
    check_backbone(tree_of("((a:1,b:1):1,c:1);"), seqs, set())


def test_seeds_needed_replaces_only_the_tail():
    hits = [1, 0, 0, 1, 0, 0, 0, 0]
    assert seeds_needed(hits, n_available=10, target=2) == 0
    assert seeds_needed(hits, n_available=10, target=4) == 2
    assert seeds_needed([0, 0, 0, 1], n_available=10, target=2) == 2


def test_seeds_needed_loses_hits_in_the_replaced_tail():
    assert seeds_needed([1, 1, 1, 0], n_available=10, target=4) == 1
    assert seeds_needed([0, 1, 1, 1], n_available=10, target=4) == 4


def test_seeds_needed_raises_when_short_of_candidates():
    with pytest.raises(ValueError, match="seed"):
        seeds_needed([0, 0, 0], n_available=1, target=3)


def test_bin_rows_cover_every_bin_and_flag_small_ones():
    rows = bin_rows("B10k", "V4-150", ["99-99.9"] * 150 + [TOP_BIN] * 3)
    frame = pd.DataFrame(rows)
    assert list(frame["bin"]) == [*BIN_LABELS, TOP_BIN]
    assert set(frame["kind"]) == {"bin"}
    assert set(frame["backbone"]) == {"B10k"}
    assert set(frame["region"]) == {"V4-150"}
    by_bin = frame.set_index("bin")
    assert by_bin.loc["99-99.9", ["n", "flagged"]].tolist() == [150, 0]
    assert by_bin.loc[TOP_BIN, ["n", "flagged"]].tolist() == [3, 1]
    assert by_bin.loc["<75", ["n", "flagged"]].tolist() == [0, 1]
    assert (frame["path"] == "").all()
    assert (frame["sha256"] == "").all()


def test_file_row_schema():
    row = file_row("data/gg2/Bfull/tree.nwk", "Bfull", 12, "ab" * 32)
    assert row == {
        "kind": "file",
        "backbone": "Bfull",
        "region": "",
        "bin": "",
        "path": "data/gg2/Bfull/tree.nwk",
        "n": 12,
        "sha256": "ab" * 32,
        "flagged": "",
    }
    assert list(row) == list(MANIFEST_COLUMNS)
    assert MANIFEST_COLUMNS == (
        "kind", "backbone", "region", "bin", "path", "n", "sha256", "flagged"
    )  # fmt: skip


def test_region_bins_on_fixture_counts_each_query_once():
    queries = dict(read_fasta(FIXTURE / "queries.aligned.fasta"))
    backbone = dict(read_fasta(FIXTURE / "aligned-dna-sequences.fasta"))
    ungapped = {name: ungap(row) for name, row in queries.items()}
    bins = region_bins(ungapped, queries, backbone, "full-length")
    assert len(bins) == len(queries)
    assert set(bins) <= {*BIN_LABELS, TOP_BIN}
    v4 = region_bins(ungapped, queries, backbone, "V4-150")
    assert 0 < len(v4) <= len(queries)


def lineages_by_cut(tree, caps):
    """Rank r taxon = the highest ancestor of each leaf with at most caps[r] leaves."""
    size = {}
    for node in tree.traverse_postorder():
        size[node] = 1 if node.is_leaf() else sum(size[c] for c in node.children)
    lineages = {}
    for leaf in tree.traverse_leaves():
        taxa = []
        for rank, cap in enumerate(caps):
            top = leaf
            while top.parent is not None and size[top.parent] <= cap:
                top = top.parent
            taxa.append(f"r{rank}_{id(top)}")
        lineages[leaf.label] = tuple(taxa)
    return lineages


def small_world(seed=3, n_leaves=200, width=60):
    rng = np.random.default_rng(seed)
    tree = coalescent_tree(n_leaves, rng, zero_fraction=0.0)
    bases = np.array(list("ACGT"))
    seqs = {
        leaf.label: "".join(bases[rng.integers(4, size=width)])
        for leaf in tree.traverse_leaves()
    }
    return tree, seqs, lineages_by_cut(tree, (60, 25, 10))


def draw_groups(tree, lineages, seed=1):
    rng = np.random.default_rng(seed)
    groups = choose_groups(
        tree, lineages, {0: 1, 1: 1, 2: 2}, 8, max_leaves=60, rng=rng, min_leaves=3
    )
    return groups, sample_queries(groups, rng, per_clade=5)


def test_build_backbones_acceptance():
    tree, seqs, lineages = small_world()
    groups, queries = draw_groups(tree, lineages)
    held = {m for g in groups for m in g.members}
    kept = [leaf.label for leaf in tree.traverse_leaves() if leaf.label not in held]
    seqs[kept[40]] = seqs[kept[3]]

    result = build_backbones(tree, seqs, held, b10k_size=50, seed=7)

    full = {leaf.label for leaf in result.full_tree.traverse_leaves()}
    small = {leaf.label for leaf in result.small_tree.traverse_leaves()}
    assert result.duplicates == {kept[40]}
    assert full == set(kept) - result.duplicates
    assert full.isdisjoint(held)
    assert small.isdisjoint(held)
    assert small <= full
    assert len(small) == 50
    assert set(result.full_seqs) == full
    assert set(result.small_seqs) == small
    check_backbone(result.full_tree, result.full_seqs, held)
    check_backbone(result.small_tree, result.small_seqs, held)
    assert set(queries) <= held
    assert all(q in seqs for q in queries)

    for backbone in (result.full_tree, result.small_tree):
        keep = {leaf.label for leaf in backbone.traverse_leaves()} | set(queries)
        reference = reference_tree(tree, keep)
        true_positions(reference.newick(), queries, prepare_backbone(backbone.newick()))

    rows = heldout_rows(groups, queries, seeds=set())
    assert {r["leaf"] for r in rows} == held
    assert len(rows) == len(held)
    assert {r["leaf"] for r in rows if r["query"] == 1} == set(queries)
    assert list(rows[0]) == [
        "leaf", "group", "group_type", "rank", "taxon", "query", "seeded"
    ]  # fmt: skip


def test_build_backbones_is_deterministic_in_seed():
    tree, seqs, lineages = small_world()
    groups, _ = draw_groups(tree, lineages)
    held = {m for g in groups for m in g.members}
    first = build_backbones(tree, seqs, held, b10k_size=40, seed=5)
    again = build_backbones(tree, seqs, held, b10k_size=40, seed=5)
    other = build_backbones(tree, seqs, held, b10k_size=40, seed=6)
    assert list(first.small_seqs) == list(again.small_seqs)
    assert list(first.small_seqs) != list(other.small_seqs)


def test_heldout_rows_mark_seeds():
    tree, seqs, lineages = small_world()
    groups, queries = draw_groups(tree, lineages)
    singleton = next(g.members[0] for g in groups if g.group_type == "singleton")
    rows = heldout_rows(groups, queries, seeds={singleton})
    assert {r["leaf"] for r in rows if r["seeded"] == 1} == {singleton}


def test_choose_seeded_groups_replaces_tail_with_near_identical_leaves():
    rng = np.random.default_rng(11)
    tree = coalescent_tree(400, rng, zero_fraction=0.0)
    bases = np.array(list("ACGT"))
    seqs = {
        leaf.label: "".join(bases[rng.integers(4, size=200)])
        for leaf in tree.traverse_leaves()
    }
    labels = list(seqs)
    for i in range(0, 100, 5):
        seqs[labels[200 + i]] = seqs[labels[i]]
    lineages = lineages_by_cut(tree, (120, 50, 20))
    options = {
        "rank_counts": {0: 1, 1: 1},
        "n_queries": 30,
        "max_leaves": 120,
        "min_leaves": 5,
        "per_clade": 5,
        "seed": 2,
        "n_candidates": 100,
    }

    def count_hits(groups):
        held = {m for g in groups for m in g.members}
        singles = {
            g.members[0]: seqs[g.members[0]]
            for g in groups
            if g.group_type == "singleton"
        }
        return sum(near_identical(singles, prune_alignment(seqs, held)))

    unseeded, _, none = choose_seeded_groups(tree, seqs, lineages, target=0, **options)
    assert none == set()
    base = count_hits(unseeded)
    target = base + 3
    groups, queries, seeds = choose_seeded_groups(
        tree, seqs, lineages, target=target, **options
    )
    assert len(queries) == 30
    assert len(seeds) >= 1
    assert count_hits(groups) > base
    assert [g for g in groups if g.group_type == "clade"] == [
        g for g in unseeded if g.group_type == "clade"
    ]


def make_qza(path, members, root="u-1"):
    with zipfile.ZipFile(path, "w") as z:
        for name, data in members.items():
            z.writestr(f"{root}/{name}", data)


def run_unzip(qza, out, *members):
    return subprocess.run(
        [sys.executable, str(SCRIPTS / "unzip_qza.py"), str(qza), str(out), *members],
        capture_output=True,
        text=True,
    )


def test_unzip_qza_extracts_requested_members(tmp_path):
    qza = tmp_path / "a.qza"
    make_qza(qza, {"data/tree.nwk": "(a,b);", "data/other.txt": "x", "VERSION": "5"})
    out = tmp_path / "out"
    result = run_unzip(qza, out, "tree.nwk")
    assert result.returncode == 0, result.stderr
    assert (out / "tree.nwk").read_text() == "(a,b);"
    assert not (out / "other.txt").exists()


def test_unzip_qza_fails_on_missing_member(tmp_path):
    qza = tmp_path / "a.qza"
    make_qza(qza, {"data/tree.nwk": "(a,b);"})
    assert run_unzip(qza, tmp_path / "out", "missing.txt").returncode != 0


@pytest.mark.parametrize("name", ["../x", "/abs"])
def test_unzip_qza_rejects_escaping_member_names(tmp_path, name):
    qza = tmp_path / "a.qza"
    make_qza(qza, {"data/tree.nwk": "(a,b);"})
    result = run_unzip(qza, tmp_path / "out", name)
    assert result.returncode != 0
    assert not (tmp_path / "x").exists()


def test_unzip_qza_rejects_two_roots(tmp_path):
    qza = tmp_path / "a.qza"
    make_qza(qza, {"data/tree.nwk": "(a,b);"})
    with zipfile.ZipFile(qza, "a") as z:
        z.writestr("u-2/data/tree.nwk", "(c,d);")
    assert run_unzip(qza, tmp_path / "out", "tree.nwk").returncode != 0


@pytest.mark.skipif(shutil.which("curl") is None, reason="curl missing")
def test_fetch_checks_sha256_before_moving_into_place(tmp_path):
    source = tmp_path / "src.txt"
    source.write_text("payload\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    out = tmp_path / "got.txt"
    fetch = SCRIPTS / "fetch.sh"
    ok = subprocess.run(
        ["bash", str(fetch), source.as_uri(), digest, str(out)], capture_output=True
    )
    assert ok.returncode == 0
    assert out.read_text() == "payload\n"
    assert not Path(f"{out}.tmp").exists()

    bad_out = tmp_path / "bad.txt"
    bad = subprocess.run(
        ["bash", str(fetch), source.as_uri(), "0" * 64, str(bad_out)],
        capture_output=True,
    )
    assert bad.returncode != 0
    assert not bad_out.exists()
    assert not Path(f"{bad_out}.tmp").exists()
