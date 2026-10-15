import re
import runpy
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest
import treeswift

from phylobench.simulate import (
    FREQS,
    RATES,
    choose_heldout,
    gtr_matrix,
    random_tree,
    simulate_alignment,
    transition_matrices,
)

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "make_fixture.py"
FIXTURE = ROOT / "tests" / "fixtures" / "tiny"

IUPAC = {
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}  # fmt: skip
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def rc(primer):
    return primer.translate(COMPLEMENT)[::-1]


def pattern(primer):
    return re.compile("".join(IUPAC.get(c, c) for c in primer))


def read_fasta(path):
    seqs = {}
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line
    return seqs


def read_tree(name):
    return treeswift.read_tree_newick((FIXTURE / name).read_text())


def leaf_labels(node):
    return {leaf.label for leaf in node.traverse_preorder(internal=False)}


def read_truth():
    rows = [
        line.split("\t") for line in (FIXTURE / "truth.tsv").read_text().splitlines()
    ]
    assert rows[0] == ["query", "group", "group_type"]
    return rows[1:]


def bipartitions(tree):
    all_leaves = leaf_labels(tree.root)
    out = set()
    for node in tree.traverse_postorder(leaves=True, internal=True):
        if node.is_root():
            continue
        side = frozenset(leaf_labels(node))
        if "t1" in side:
            side = frozenset(all_leaves - side)
        out.add(side)
    return out


def test_regeneration_is_byte_identical(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["make_fixture.py", str(tmp_path)])
    runpy.run_path(str(SCRIPT), run_name="__main__")
    committed = {p.name for p in FIXTURE.iterdir()}
    assert {p.name for p in tmp_path.iterdir()} == committed
    for name in committed:
        assert (tmp_path / name).read_bytes() == (FIXTURE / name).read_bytes(), name


def test_total_size_under_1mb():
    assert sum(p.stat().st_size for p in FIXTURE.iterdir()) < 1_000_000


def test_names_are_consistent():
    reference = leaf_labels(read_tree("reference_tree.nwk").root)
    backbone = leaf_labels(read_tree("tree.nwk").root)
    truth = read_truth()
    queries = [row[0] for row in truth]
    backbone_fasta = read_fasta(FIXTURE / "aligned-dna-sequences.fasta")
    assert not backbone & set(queries)
    assert backbone | set(queries) == reference
    assert set(backbone_fasta) == backbone
    assert list(read_fasta(FIXTURE / "queries.fasta")) == queries
    assert list(read_fasta(FIXTURE / "queries.aligned.fasta")) == queries


FORWARD_REVERSE = [
    ("AGAGTTTGATCCTGGCTCAG", rc("GCTGCCTCCCGTAGGAGT")),
    ("CCTACGGGNBGCASCAG", rc("GACTACHVGGGTATCTAATCC")),
    ("GTGYCAGCMGCCGCGGTAA", rc("GGACTACNVGGGTWTCTAAT")),
]


@pytest.mark.parametrize("fasta", ["aligned-dna-sequences.fasta", "queries.fasta"])
def test_primer_sites_match_once_in_order(fasta):
    for name, seq in read_fasta(FIXTURE / fasta).items():
        for forward, reverse in FORWARD_REVERSE:
            f_hits = list(pattern(forward).finditer(seq))
            r_hits = list(pattern(reverse).finditer(seq))
            assert len(f_hits) == 1 and len(r_hits) == 1, (name, forward, reverse)
            assert f_hits[0].start() < r_hits[0].start()


def test_raxml_info_layout_and_values():
    text = (FIXTURE / "raxml-info.txt").read_text()
    assert re.search(r"^Partition: 0$", text, re.M)
    assert re.search(r"^Substitution Matrix: GTR$", text, re.M)
    alpha = re.search(
        r"^alpha\[0\]: ([\d.]+) rates\[0\] ac ag at cg ct gt: (.+)$", text, re.M
    )
    assert float(alpha[1]) == 1.0
    rates = [float(x) for x in alpha[2].split()]
    assert rates == pytest.approx(RATES, abs=1e-9) and rates[-1] == 1.0
    freqs = [
        float(x) for x in re.search(r"^Base frequencies: (.+)$", text, re.M)[1].split()
    ]
    assert freqs == pytest.approx(FREQS, abs=1e-9)
    assert sum(freqs) == pytest.approx(1.0, abs=1e-5)


def test_backbone_topology_and_distances():
    reference = read_tree("reference_tree.nwk")
    backbone = read_tree("tree.nwk")
    assert backbone.root.num_children() == 3
    assert all(n.num_children() != 1 for n in backbone.traverse_preorder())
    ref_dist = reference.distance_matrix(leaf_labels=True)
    bb_dist = backbone.distance_matrix(leaf_labels=True)
    for a, b in combinations(bb_dist, 2):
        assert bb_dist[a][b] == pytest.approx(ref_dist[a][b], abs=1e-5)


def test_clades_are_monophyletic_in_reference():
    reference = read_tree("reference_tree.nwk")
    leaves = {leaf.label: leaf for leaf in reference.traverse_leaves()}
    clades = {}
    singletons = []
    for query, group, group_type in read_truth():
        if group_type == "clade":
            clades.setdefault(group, []).append(query)
        else:
            singletons.append(query)
    assert len(clades) == 2 and len(singletons) == 4
    assert sum(len(m) for m in clades.values()) + len(singletons) == 10
    for members in clades.values():
        assert len(members) == 3
        nodes = [leaves[m] for m in members]
        common = set.intersection(
            *({id(a) for a in n.traverse_ancestors()} for n in nodes)
        )
        mrca = min(
            (a for n in nodes[:1] for a in n.traverse_ancestors() if id(a) in common),
            key=lambda a: len(leaf_labels(a)),
        )
        assert leaf_labels(mrca) == set(members)


def test_every_group_is_a_distinct_backbone_bipartition():
    reference = read_tree("reference_tree.nwk")
    backbone = read_tree("tree.nwk")
    backbone_leaves = leaf_labels(backbone.root)
    backbone_splits = bipartitions(backbone)
    leaves = {leaf.label: leaf for leaf in reference.traverse_leaves()}
    members_by_group = {}
    for query, group, _ in read_truth():
        members_by_group.setdefault(group, []).append(query)
    seen = set()
    for members in members_by_group.values():
        nodes = [leaves[m] for m in members]
        stem = nodes[0]
        while leaf_labels(stem) != set(members):
            stem = stem.parent
        assert not stem.parent.is_root()
        side = frozenset(leaf_labels(stem.parent) & backbone_leaves)
        assert side and backbone_leaves - side
        if "t1" in side:
            side = frozenset(backbone_leaves - side)
        assert side in backbone_splits
        assert side not in seen
        seen.add(side)


def test_choose_heldout_raises_when_impossible():
    rng = np.random.default_rng(0)
    tree = random_tree(110, rng)
    with pytest.raises(ValueError):
        choose_heldout(tree, (200,), 0, rng)


def test_sequences_are_clean_dna():
    backbone = read_fasta(FIXTURE / "aligned-dna-sequences.fasta")
    queries = read_fasta(FIXTURE / "queries.fasta")
    aligned = read_fasta(FIXTURE / "queries.aligned.fasta")
    for seq in [*backbone.values(), *queries.values()]:
        assert len(seq) == 1500 and set(seq) <= set("ACGT")
    assert {k: v.replace("-", "") for k, v in aligned.items()} == queries


def test_gtr_matrix_properties():
    q = np.array(gtr_matrix(RATES, FREQS))
    pi = np.array(FREQS)
    assert q.sum(axis=1) == pytest.approx(0, abs=1e-12)
    assert pi @ q == pytest.approx(0, abs=1e-6)
    assert -(pi * np.diag(q)).sum() == pytest.approx(1.0, abs=1e-12)


def test_transition_matrices_are_stochastic():
    q = gtr_matrix(RATES, FREQS)
    lengths = np.array([0.0, 0.005, 0.03, 0.055])
    p = transition_matrices(lengths, q)
    assert p.shape == (4, 4, 4, 4)
    assert p.sum(axis=-1) == pytest.approx(1.0, abs=1e-12)
    assert p.min() >= -1e-12
    assert p[0] == pytest.approx(np.broadcast_to(np.eye(4), (4, 4, 4)), abs=1e-15)


def test_scales_to_many_leaves():
    rng = np.random.default_rng(1)
    tree = random_tree(10_000, rng)
    leaves = list(tree.traverse_leaves())
    assert len({leaf.label for leaf in leaves}) == 10_000
    for node in tree.traverse_preorder(internal=True, leaves=True):
        if not node.is_root():
            units = node.edge_length * 1_000_000
            assert units == pytest.approx(round(units), abs=1e-6)
    small = simulate_alignment(random_tree(1_000, rng), 50, rng)
    assert len(small) == 1_000 and {len(s) for s in small.values()} == {50}
