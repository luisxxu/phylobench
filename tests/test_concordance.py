import math

import pytest
from tree_helpers import edge_numbers, place

from phylobench.concordance import score_concordance
from phylobench.jplace import numbered_doc

BINARY = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
ZERO = "((a:1,b:1)X:0,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"


def doc_with(placed, newick=BINARY):
    doc = numbered_doc(newick)
    numbers = edge_numbers(doc)
    for query, label in placed.items():
        place(doc, query, numbers[label])
    return doc


def concordance(placed, reference, queries, newick=BINARY, reference_newick=None):
    rows = score_concordance(
        doc_with(placed, newick),
        doc_with(reference, reference_newick or newick),
        newick,
        queries,
    )
    return {row["query"]: row["edge_dist"] for row in rows}


def test_same_placement_is_zero():
    assert concordance({"q": "a"}, {"q": "a"}, ["q"]) == {"q": 0}


def test_adjacent_edges_are_one():
    assert concordance({"q": "a"}, {"q": "b"}, ["q"]) == {"q": 1}


def test_distance_is_symmetric():
    forward = concordance({"q": "a"}, {"q": "c"}, ["q"])
    backward = concordance({"q": "c"}, {"q": "a"}, ["q"])
    assert forward == backward and forward["q"] > 1


def test_zero_length_edge_collapses_into_its_vertex():
    # X has zero length, so a and b are both adjacent to the vertex X collapses into
    assert concordance({"q": "X"}, {"q": "a"}, ["q"], ZERO) == {"q": 0}


def test_query_unplaced_in_either_doc_is_nan():
    out = concordance({"p": "a"}, {"p": "a", "q": "b"}, ["p", "q"])
    assert out["p"] == 0 and math.isnan(out["q"])
    out = concordance({"p": "a", "q": "b"}, {"p": "a"}, ["p", "q"])
    assert out["p"] == 0 and math.isnan(out["q"])


def test_both_unplaced_is_nan_without_error():
    assert math.isnan(concordance({}, {}, ["q"])["q"])


def test_rows_follow_query_order():
    rows = score_concordance(
        doc_with({"p": "a", "q": "a"}),
        doc_with({"p": "a", "q": "a"}),
        BINARY,
        ["q", "p"],
    )
    assert [row["query"] for row in rows] == ["q", "p"]


def test_unknown_query_raises():
    with pytest.raises(ValueError, match="not requested"):
        concordance({"q": "a", "x": "a"}, {"q": "a"}, ["q"])
    with pytest.raises(ValueError, match="not requested"):
        concordance({"q": "a"}, {"q": "a", "x": "a"}, ["q"])


def test_duplicate_queries_raise():
    with pytest.raises(ValueError, match="duplicate"):
        concordance({"q": "a"}, {"q": "a"}, ["q", "q"])


def test_mismatched_leaves_raise():
    other = "((a:1,b:1)X:1,((c:1,d:1)Y:1,g:1)Z:1,f:1)R;"
    with pytest.raises(ValueError, match="leaves differ"):
        score_concordance(
            doc_with({"q": "a"}), doc_with({"q": "a"}, other), BINARY, ["q"]
        )
