from itertools import count
from pathlib import Path

import numpy as np
import pytest
import treeswift

from phylobench.accuracy import prepare_backbone, true_positions
from phylobench.heldout import (
    HeldOutGroup,
    choose_groups,
    monophyletic_taxa,
    prune_alignment,
    prune_tree,
    sample_queries,
)
from phylobench.simulate import format_newick

FIXTURE = Path(__file__).parent / "fixtures" / "tiny"


def read_tree(newick):
    return treeswift.read_tree_newick(newick)


def leaf_set(node):
    return {leaf.label for leaf in node.traverse_leaves()}


def balanced(labels):
    if len(labels) == 1:
        return f"{labels[0]}:1"
    mid = len(labels) // 2
    return f"({balanced(labels[:mid])},{balanced(labels[mid:])}):1"


def taxonomy_tree():
    """root -> 4 halves -> 2 orders each -> 4 families -> 3 genera of 30 leaves."""
    serial = count()
    lineages = {}

    def genus(order, family, g):
        labels = [f"l{next(serial)}" for _ in range(30)]
        for label in labels:
            lineages[label] = (
                f"O{order}",
                f"F{order}.{family}",
                f"G{order}.{family}.{g}",
            )
        return balanced(labels)

    def family(order, f):
        return "(" + ",".join(genus(order, f, g) for g in range(3)) + "):1"

    def order(o):
        return "(" + ",".join(family(o, f) for f in range(4)) + "):1"

    halves = [f"({order(2 * h)},{order(2 * h + 1)}):1" for h in range(4)]
    return read_tree(f"({','.join(halves)});"), lineages


def parent_leaf_set(tree, group):
    return leaf_set(tree.mrca(set(group.members)).parent)


RANKS = {0: 1, 1: 2, 2: 2}


def pick(tree, lineages, seed=0, n_singletons=10, rank_counts=None, **kwargs):
    return choose_groups(
        tree,
        lineages,
        rank_counts or RANKS,
        n_singletons,
        max_leaves=400,
        rng=np.random.default_rng(seed),
        **kwargs,
    )


def test_prune_tree_contracts_degree_two_root():
    tree = read_tree("((a:1,b:1,x:1):2,(c:1,d:1):3);")
    out = prune_tree(tree, {"x"})
    assert format_newick(out) == (
        "((c:1.000000,d:1.000000):5.000000,a:1.000000,b:1.000000);\n"
    )


def test_prune_tree_clears_stale_root_edge_and_adds_unary_lengths():
    tree = read_tree("(((a:1,b:1):2,x:1):3,c:1,d:1):7;")
    out = prune_tree(tree, {"x"})
    assert out.root.edge_length is None
    text = format_newick(out)
    assert text == "(c:1.000000,d:1.000000,(a:1.000000,b:1.000000):5.000000);\n"
    assert "[" not in text


def test_prune_tree_does_not_mutate_input_and_rejects_unknown_leaf():
    tree = read_tree("((a:1,b:1):1,c:1,d:1);")
    before = format_newick(tree)
    prune_tree(tree, {"a"})
    assert format_newick(tree) == before
    with pytest.raises(ValueError, match="nope"):
        prune_tree(tree, {"nope"})


def read_truth():
    rows = [
        line.split("\t") for line in (FIXTURE / "truth.tsv").read_text().splitlines()
    ]
    return [(r[0], r[1], r[2]) for r in rows[1:]]


def test_fixture_prune_reproduces_committed_backbone():
    reference_text = (FIXTURE / "reference_tree.nwk").read_text()
    queries = {q for q, _, _ in read_truth()}
    out = prune_tree(read_tree(reference_text), queries)
    text = format_newick(out)
    assert text == (FIXTURE / "tree.nwk").read_text()
    assert queries.isdisjoint(leaf_set(out.root))
    assert all(len(n.children) != 1 for n in out.traverse_preorder())


def test_fixture_every_group_has_a_backbone_stem():
    reference_text = (FIXTURE / "reference_tree.nwk").read_text()
    truth = read_truth()
    queries = [q for q, _, _ in truth]
    backbone = prepare_backbone((FIXTURE / "tree.nwk").read_text())
    node, _ = true_positions(reference_text, queries, backbone)
    by_group = {}
    for (_, group, _), n in zip(truth, node.tolist(), strict=True):
        by_group.setdefault(group, set()).add(n)
    assert all(len(nodes) == 1 for nodes in by_group.values())
    assert len({next(iter(n)) for n in by_group.values()}) == len(by_group)


def test_monophyletic_taxa_excludes_split_taxa():
    tree = read_tree("(((a,b),c),(d,e));")
    lineages = {
        "a": ("p1", "x"),
        "b": ("p1", "x"),
        "c": ("p2", "y"),
        "d": ("p1", "z"),
        "e": ("p3", "w"),
    }
    top = monophyletic_taxa(tree, lineages, 0)
    assert set(top) == {"p2", "p3"}
    assert {leaf_set(n).pop() for n in top.values()} == {"c", "e"}
    second = monophyletic_taxa(tree, lineages, 1)
    assert leaf_set(second["x"]) == {"a", "b"}


def test_monophyletic_taxa_leaf_without_lineage_breaks_monophyly():
    tree = read_tree("(((a,b,c),d),e);")
    lineages = {"a": ("t",), "b": ("t",), "d": ("u",), "e": ("v",)}
    assert set(monophyletic_taxa(tree, lineages, 0)) == {"u", "v"}


def test_monophyletic_taxa_ignores_leaves_shorter_than_rank():
    tree = read_tree("((a,b),c);")
    lineages = {"a": ("t", "x"), "b": ("t",), "c": ("u", "y")}
    assert set(monophyletic_taxa(tree, lineages, 1)) == {"x", "y"}


def test_choose_groups_shapes_and_disjointness():
    tree, lineages = taxonomy_tree()
    groups = pick(tree, lineages)
    clades = [g for g in groups if g.group_type == "clade"]
    singles = [g for g in groups if g.group_type == "singleton"]
    assert [g.rank for g in clades] == [0, 1, 1, 2, 2]
    assert len(singles) == 10
    assert all(g.rank is None and len(g.members) == 1 for g in singles)
    assert all(25 <= len(g.members) <= 400 for g in clades)
    assert len({g.group_id for g in groups}) == len(groups)
    all_members = [m for g in groups for m in g.members]
    assert len(all_members) == len(set(all_members))
    for g in clades:
        taxon_leaves = {
            leaf for leaf, lin in lineages.items() if lin[g.rank] == g.taxon
        }
        assert set(g.members) == taxon_leaves


def test_choose_groups_stem_parents_untouched_and_not_root():
    tree, lineages = taxonomy_tree()
    groups = pick(tree, lineages)
    for g in groups:
        stem_parent = tree.mrca(set(g.members)).parent
        assert not stem_parent.is_root()
        others = {m for o in groups if o is not g for m in o.members}
        assert leaf_set(stem_parent).isdisjoint(others)


def test_choose_groups_raises_when_exhausted():
    tree, lineages = taxonomy_tree()
    with pytest.raises(ValueError):
        pick(tree, lineages, rank_counts={2: 100})
    with pytest.raises(ValueError):
        pick(tree, lineages, n_singletons=100_000)
    with pytest.raises(ValueError):
        choose_groups(
            tree, lineages, {0: 1}, 0, max_leaves=100, rng=np.random.default_rng(0)
        )
    with pytest.raises(ValueError):
        choose_groups(
            tree,
            lineages,
            {2: 1},
            0,
            max_leaves=400,
            rng=np.random.default_rng(0),
            min_leaves=31,
        )


def test_choose_groups_is_deterministic_per_seed():
    tree, lineages = taxonomy_tree()
    assert pick(tree, lineages, seed=5) == pick(tree, lineages, seed=5)
    assert pick(tree, lineages, seed=5) != pick(tree, lineages, seed=6)


def test_sample_queries_caps_clades_and_keeps_all_singletons():
    groups = [
        HeldOutGroup("clade1", "clade", 2, "G", tuple(f"c{i}" for i in range(30))),
        HeldOutGroup("clade2", "clade", 2, "H", ("d0", "d1", "d2")),
        HeldOutGroup("singleton1", "singleton", None, "s0", ("s0",)),
        HeldOutGroup("singleton2", "singleton", None, "s1", ("s1",)),
    ]
    picked = sample_queries(groups, np.random.default_rng(0), 25)
    assert len(picked) == len(set(picked)) == 25 + 3 + 2
    assert sum(q.startswith("c") for q in picked) == 25
    assert {"d0", "d1", "d2", "s0", "s1"} <= set(picked)
    assert picked == sample_queries(groups, np.random.default_rng(0), 25)


def test_prune_alignment_keeps_columns_and_rejects_missing():
    seqs = {"a": "AC-T", "b": "ACGT", "c": "A--T"}
    assert prune_alignment(seqs, {"b"}) == {"a": "AC-T", "c": "A--T"}
    with pytest.raises(ValueError, match="zz"):
        prune_alignment(seqs, {"zz"})
