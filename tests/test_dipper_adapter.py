import math

import pytest
from tree_helpers import edge_numbers

from phylobench.accuracy import score_jplace
from phylobench.dipper_adapter import dipper_to_jplace
from phylobench.jplace import REQUIRED_FIELDS, edge_numbered_tree, top_placements

BACKBONE = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"


def convert(dipper, queries, backbone=BACKBONE):
    return dipper_to_jplace(dipper, backbone, queries)


def placed_labels(doc):
    numbers = edge_numbers(doc)
    by_edge = {n: label for label, n in numbers.items()}
    return {p.query: by_edge[p.edge_num] for p in top_placements(doc)}


def scores(doc, dipper, queries, backbone=BACKBONE):
    rows = score_jplace(doc, backbone, dipper, queries)
    return {r["query"]: r for r in rows}


def test_single_query_on_an_edge():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,(e:0.5,q:0.1):0.5)Z:1,f:1)R;"
    doc = convert(dipper, ["q"])
    assert placed_labels(doc) == {"q": "e"}
    assert scores(doc, dipper, ["q"])["q"]["edge_dist"] == 0


def test_query_cherry_maps_to_the_edge_where_the_clade_attaches():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,(f:0.5,(q1:0.1,q2:0.1):0.2):0.5)R;"
    doc = convert(dipper, ["q1", "q2"])
    assert placed_labels(doc) == {"q1": "f", "q2": "f"}


def test_two_queries_on_the_same_backbone_edge():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,((e:0.3,q1:0.1):0.2,q2:0.1):0.5)Z:1,f:1)R;"
    doc = convert(dipper, ["q1", "q2"])
    assert placed_labels(doc) == {"q1": "e", "q2": "e"}
    assert all(r["edge_dist"] == 0 for r in scores(doc, dipper, ["q1", "q2"]).values())


def test_query_attached_to_a_query():
    chain = (
        "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,"
        "(f:0.5,(q1:0.1,(q2:0.1,x:0.1):0.1):0.2):0.5)R;"
    )
    doc = convert(chain, ["q1", "q2", "x"])
    assert placed_labels(doc) == {"q1": "f", "q2": "f", "x": "f"}


def test_rerooted_at_a_backbone_leaf():
    rerooted = "(f:0.5,q:0.1,(((c:1,d:1)Y:1,e:1)Z:1,(a:1,b:1)X:1):0.5)R;"
    doc = convert(rerooted, ["q"])
    assert placed_labels(doc) == {"q": "f"}
    assert scores(doc, rerooted, ["q"])["q"]["edge_dist"] == 0


def test_internal_vertex_attachment_scores_zero_against_itself():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    dipper = dipper.replace("(c:1,d:1)Y:1", "(c:1,d:1,q:0.1)Y:1")
    doc = convert(dipper, ["q"])
    assert scores(doc, dipper, ["q"])["q"]["edge_dist"] == 0
    assert placed_labels(doc) == {"q": "Y"}


def test_root_vertex_attachment_scores_zero_against_itself():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1,q:0.1)R;"
    doc = convert(dipper, ["q"])
    assert scores(doc, dipper, ["q"])["q"]["edge_dist"] == 0


def test_vertex_collapse_pessimism_up_to_one():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;".replace(
        "(c:1,d:1)Y:1", "(c:1,d:1,q:0.1)Y:1"
    )
    doc = convert(dipper, ["q"])
    truth_on_child_edge = "((a:1,b:1)X:1,(((c:1,q:0.1):0.5,d:1)Y:1,e:1)Z:1,f:1)R;"
    row = scores(doc, truth_on_child_edge, ["q"])["q"]
    assert row["edge_dist"] == 1


def test_zero_length_backbone_edge():
    backbone = "((a:1,b:1)X:0,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    dipper = "((a:1,b:1)X:0,((c:1,d:1)Y:1,(e:0.5,q:0.1):0.5)Z:1,f:1)R;"
    doc = convert(dipper, ["q"], backbone)
    assert scores(doc, dipper, ["q"], backbone)["q"]["edge_dist"] == 0


def test_bifurcating_root_backbone_uses_the_lowest_edge_number():
    backbone = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1)R;"
    dipper = "(((a:1,b:1)X:0.5,q:0.1):0.5,((c:1,d:1)Y:1,e:1)Z:1)R;"
    doc = convert(dipper, ["q"], backbone)
    numbers = edge_numbers(doc)
    assert doc["placements"][0]["p"][0][0] == numbers["X"] < numbers["Z"]
    assert scores(doc, dipper, ["q"], backbone)["q"]["edge_dist"] == 0


def test_changed_topology_raises():
    dipper = "((a:1,c:1)X:1,((b:1,d:1)Y:1,e:1)Z:1,(f:1,q:0.1):0)R;"
    with pytest.raises(ValueError, match="differs"):
        convert(dipper, ["q"])


def test_missing_backbone_leaf_raises():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,q:0.1)Z:1,f:1)R;"
    with pytest.raises(ValueError):
        convert(dipper, ["q"])


def test_extra_leaf_raises():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,(f:1,q:0.1,stray:0.1):0)R;"
    with pytest.raises(ValueError):
        convert(dipper, ["q"])


def test_query_named_like_a_backbone_leaf_raises():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    with pytest.raises(ValueError, match="backbone"):
        convert(dipper, ["a"])


def test_duplicate_queries_raise():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,(f:1,q:0.1):0)R;"
    with pytest.raises(ValueError, match="duplicate"):
        convert(dipper, ["q", "q"])


def test_absent_query_scores_nan_and_others_are_scored():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,(e:0.5,q:0.1):0.5)Z:1,f:1)R;"
    reference = "((a:1,b:1)X:1,((c:1,d:1)Y:1,(e:0.5,q:0.1,gone:0.1):0.5)Z:1,f:1)R;"
    doc = convert(dipper, ["q", "gone"])
    rows = scores(doc, reference, ["q", "gone"])
    assert rows["q"]["edge_dist"] == 0
    assert math.isnan(rows["gone"]["edge_dist"])


def test_output_shape():
    dipper = "((a:1,b:1)X:1,((c:1,d:1)Y:1,(e:0.5,q:0.1):0.5)Z:1,f:1)R;"
    doc = convert(dipper, ["q"])
    assert doc["version"] == 3
    assert set(REQUIRED_FIELDS) <= set(doc["fields"])
    row = doc["placements"][0]["p"][0]
    assert row[1] == 1.0
    assert math.isnan(row[2]) and math.isnan(row[3])
    tree, _ = edge_numbered_tree(doc["tree"])
    assert {leaf.label for leaf in tree.traverse_leaves()} == set("abcdef")
    assert [p.query for p in top_placements(doc)] == ["q"]
    assert "q" not in doc["tree"]
