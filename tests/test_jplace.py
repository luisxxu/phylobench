import pytest

from phylobench.jplace import Placement, edge_numbered_tree, top_placements

PPLACER_FIELDS = ["distal_length", "edge_num", "like_weight_ratio", "pendant_length"]
EPA_FIELDS = [
    "edge_num",
    "likelihood",
    "like_weight_ratio",
    "distal_length",
    "pendant_length",
]


def doc(fields, rows, **record):
    return {
        "version": 3,
        "fields": fields,
        "tree": "(a,b);",
        "placements": [{"p": rows, **record}],
    }


def test_field_order_does_not_matter():
    pplacer = doc(PPLACER_FIELDS, [[0.2, 7, 0.9, 0.3]], n=["q"])
    epa = doc(EPA_FIELDS, [[7, -10.0, 0.9, 0.2, 0.3]], n=["q"])
    expected = [Placement("q", 7, 0.9, 0.2, 0.3)]
    assert top_placements(pplacer) == top_placements(epa) == expected


def test_top_lwr_is_not_necessarily_first():
    rows = [[0.0, 1, 0.1, 0.0], [0.0, 2, 0.7, 0.0], [0.0, 3, 0.2, 0.0]]
    assert top_placements(doc(PPLACER_FIELDS, rows, n=["q"]))[0].edge_num == 2


def test_lwr_tie_keeps_first_row():
    rows = [[0.0, 5, 0.5, 0.0], [0.0, 6, 0.5, 0.0]]
    assert top_placements(doc(PPLACER_FIELDS, rows, n=["q"]))[0].edge_num == 5


def test_multi_name_record_gives_one_row_per_name():
    rows = [[0.1, 4, 1.0, 0.2]]
    result = top_placements(doc(PPLACER_FIELDS, rows, n=["a", "b"]))
    assert [(p.query, p.edge_num) for p in result] == [("a", 4), ("b", 4)]


def test_nm_records():
    rows = [[0.1, 4, 1.0, 0.2]]
    result = top_placements(doc(PPLACER_FIELDS, rows, nm=[["a", 3], ["b", 1]]))
    assert [p.query for p in result] == ["a", "b"]


def test_empty_placement_list_is_skipped():
    assert top_placements(doc(PPLACER_FIELDS, [], n=["q"])) == []


def test_missing_required_field_raises():
    with pytest.raises(ValueError, match="missing"):
        top_placements(doc(["edge_num", "like_weight_ratio"], [[1, 1.0]], n=["q"]))


def test_unsupported_version_raises():
    bad = doc(PPLACER_FIELDS, [[0.1, 4, 1.0, 0.2]], n=["q"]) | {"version": 2}
    with pytest.raises(ValueError, match="version"):
        top_placements(bad)


def test_duplicate_query_raises():
    d = doc(PPLACER_FIELDS, [[0.1, 4, 1.0, 0.2]], n=["q"])
    d["placements"].append({"p": [[0.1, 5, 1.0, 0.2]], "n": ["q"]})
    with pytest.raises(ValueError, match="duplicate"):
        top_placements(d)


def labels_and_numbers(newick):
    tree, nodes = edge_numbered_tree(newick)
    return {number: node.label for number, node in nodes.items()}, tree


def test_internal_edge_number_after_length():
    numbers, tree = labels_and_numbers("((a:1{0},b:1{1})x:0.1{3},c:2{2});")
    assert numbers == {0: "a", 1: "b", 3: "x", 2: "c"}
    assert tree.root.label is None


def test_unlabelled_internal_node_gets_none():
    numbers, _ = labels_and_numbers("((a:1{0},b:1{1}):0.1{3},c:2{2});")
    assert numbers[3] is None


def test_root_edge_number_without_length():
    numbers, _ = labels_and_numbers("(a:1{0},b:1{1}){4};")
    assert numbers[4] is None


def test_root_with_zero_length_and_number():
    numbers, _ = labels_and_numbers("((a:1{0},b:1{1}):0{5},c:2{2}):0{6};")
    assert numbers[6] is None and numbers[5] is None


@pytest.mark.parametrize("length", ["1e-05", "-0.5", "2.0E+1"])
def test_unusual_length_formats(length):
    tree, nodes = edge_numbered_tree(f"(a:{length}{{0}},b:1{{1}});")
    assert nodes[0].label == "a"
    assert nodes[0].edge_length == float(length)


def test_quoted_label_with_space():
    _, nodes = edge_numbered_tree("('a b':1{0},c:1{1});")
    assert nodes[0].label == "a b"


def test_label_containing_at_sign():
    _, nodes = edge_numbered_tree("(x@y:1{0},c:1{1});")
    assert nodes[0].label == "x@y"


def test_leading_rooted_comment():
    _, nodes = edge_numbered_tree("[&R] (a:1{0},b:1{1});")
    assert set(nodes) == {0, 1}


def test_duplicate_edge_number_raises():
    with pytest.raises(ValueError, match="duplicate"):
        edge_numbered_tree("(a:1{0},b:1{0});")


def test_unnumbered_label_with_at_sign_is_kept():
    tree, nodes = edge_numbered_tree("(a:1{0},b:1{1})foo@bar;")
    assert set(nodes) == {0, 1}
    assert tree.root.label == "foo@bar"
