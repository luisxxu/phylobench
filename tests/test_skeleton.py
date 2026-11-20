import math
from pathlib import Path

import numpy as np
import pytest
import treeswift
from tree_helpers import edge_numbers, place

from phylobench.accuracy import (
    canonical,
    edge_numbered_tree,
    prepare_backbone,
    score_jplace,
    true_positions,
)
from phylobench.jplace import numbered_doc
from phylobench.regions import read_fasta
from phylobench.skeleton import skeleton_map, skeleton_newick, skeleton_scores

TINY = Path(__file__).parent / "fixtures" / "tiny"
G = {"a", "b", "c", "d", "e"}
# n1 and n2 are not G leaves; n2 hangs off the root.
BACKBONE = "(((a:1,b:1)X:1,c:1)W:1,((d:1,n1:1)Y:1,e:1)Z:1,n2:1)R;"
REFERENCE = "(((a:1,b:1)X:1,c:1)W:1,(((d:1,q:1)Dq:1,n1:1)Y:1,e:1)Z:1,n2:1)R;"
INDUCED = "(((a:1,b:1)X:1,c:1)W:1,((d:1,q:1)Dq:1,e:1)Z:1)R;"


def node_of(backbone, *leaves):
    sig = 0
    for leaf in leaves:
        sig ^= backbone.hashes[leaf]
    return backbone.key_to_node[canonical(sig, backbone.total)]


def skeleton_of(newick, keep):
    return prepare_backbone(skeleton_newick(newick, keep))


def depths(newick):
    tree = treeswift.read_tree_newick(newick)
    return {
        leaf.label: tree.distance_between(tree.root, leaf)
        for leaf in tree.traverse_leaves()
    }


def score(backbone_newick, reference, keep, queries, placed):
    backbone = prepare_backbone(backbone_newick)
    skeleton = skeleton_of(backbone_newick, keep)
    doc = numbered_doc(backbone_newick)
    numbers = edge_numbers(doc)
    for query, label in placed.items():
        place(doc, query, numbers[label])
    true_node, true_weight = true_positions(reference, queries, skeleton)
    rows = skeleton_scores(
        doc,
        backbone,
        skeleton,
        skeleton_map(backbone_newick, skeleton),
        true_node,
        true_weight,
        queries,
    )
    return {row["query"]: row["edge_dist"] for row in rows}, doc


def test_skeleton_newick_keeps_leaves_suppresses_unary_and_sums_lengths():
    out = skeleton_newick(BACKBONE, G)
    assert set(depths(out)) == G
    assert depths(out)["d"] == depths(BACKBONE)["d"]
    assert depths(out)["a"] == depths(BACKBONE)["a"]
    tree = treeswift.read_tree_newick(out)
    assert all(len(n.children) != 1 for n in tree.traverse_preorder())
    d = next(leaf for leaf in tree.traverse_leaves() if leaf.label == "d")
    assert d.edge_length == 2.0


def test_skeleton_newick_needs_four_leaves():
    with pytest.raises(ValueError, match="at least 4"):
        skeleton_newick(BACKBONE, {"a", "b", "c"})


def test_map_sends_non_g_nodes_to_the_g_edge_they_hang_from():
    backbone = prepare_backbone(BACKBONE)
    skeleton = skeleton_of(BACKBONE, G)
    node_map = skeleton_map(BACKBONE, skeleton)
    d_edge = node_of(skeleton, "d")
    assert node_map[node_of(backbone, "n1")] == d_edge
    assert node_map[node_of(backbone, "d", "n1")] == d_edge
    assert node_map[node_of(backbone, "d")] == d_edge
    assert node_map[node_of(backbone, "a", "b")] == node_of(skeleton, "a", "b")
    z_edge = node_of(skeleton, "d", "e")
    assert node_map[node_of(backbone, "d", "e", "n1")] == z_edge
    assert node_map[node_of(backbone, "a", "b", "c")] == z_edge
    assert node_map[node_of(backbone, "n2")] == 0
    assert node_map[0] == 0


def test_placement_on_a_non_g_pendant_scores_on_the_edge_it_hangs_from():
    queries = ["q"]
    dist, doc = score(BACKBONE, REFERENCE, G, queries, {"q": "n1"})
    assert dist["q"] == 0
    full = {r["query"]: r for r in score_jplace(doc, BACKBONE, REFERENCE, queries)}
    assert full["q"]["edge_dist"] >= 1


def test_unplaced_query_is_nan():
    dist, _ = score(BACKBONE, REFERENCE, G, ["q"], {})
    assert math.isnan(dist["q"])


def test_skeleton_truth_equals_truth_on_the_hand_built_induced_tree():
    skeleton = skeleton_of(BACKBONE, G)
    from_full = true_positions(REFERENCE, ["q"], skeleton)
    from_induced = true_positions(INDUCED, ["q"], skeleton)
    assert all(
        np.array_equal(x, y) for x, y in zip(from_full, from_induced, strict=True)
    )


@pytest.mark.parametrize(
    "newick",
    [
        "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;",
        "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1)R;",
    ],
)
def test_identity_map_when_every_leaf_is_kept(newick):
    keep = {"a", "b", "c", "d", "e", "f"}
    skeleton = skeleton_of(newick, keep)
    assert skeleton_map(newick, skeleton).tolist() == list(range(len(skeleton.sig)))


def test_keeping_every_leaf_reproduces_score_jplace_on_the_tiny_fixture():
    backbone_newick = (TINY / "tree.nwk").read_text()
    reference = (TINY / "reference_tree.nwk").read_text()
    queries = [name for name, _ in read_fasta(TINY / "queries.fasta")]
    doc = numbered_doc(backbone_newick)
    n_edges = len(edge_numbered_tree(doc["tree"])[1])
    for i, query in enumerate(queries):
        place(doc, query, i * 7 % n_edges)
    keep = set(prepare_backbone(backbone_newick).hashes)
    backbone = prepare_backbone(backbone_newick)
    skeleton = skeleton_of(backbone_newick, keep)
    true_node, true_weight = true_positions(reference, queries, skeleton)
    rows = skeleton_scores(
        doc,
        backbone,
        skeleton,
        skeleton_map(backbone_newick, skeleton),
        true_node,
        true_weight,
        queries,
    )
    expected = score_jplace(doc, backbone_newick, reference, queries)
    assert [r["edge_dist"] for r in rows] == [r["edge_dist"] for r in expected]
