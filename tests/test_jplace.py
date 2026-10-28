import math

import pytest

from phylobench.accuracy import jplace_edge_nodes, prepare_backbone, score_jplace
from phylobench.jplace import (
    Placement,
    edge_numbered_tree,
    from_sepp,
    parse_jplace,
    top_placements,
)

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


def test_all_non_finite_lwr_record_is_skipped():
    rows = [[0.0, 1, math.nan, 2.0], [0.0, 2, math.nan, 2.0]]
    assert top_placements(doc(PPLACER_FIELDS, rows, n=["q"])) == []


def test_non_finite_lwr_rows_are_ignored():
    rows = [[0.0, 1, math.nan, 0.0], [0.0, 2, 0.4, 0.0], [0.0, 3, -math.inf, 0.0]]
    assert top_placements(doc(PPLACER_FIELDS, rows, n=["q"]))[0].edge_num == 2


def test_parse_jplace_reads_pplacer_non_finite_numbers():
    text = (
        '{"fields": ["edge_num", "like_weight_ratio", "likelihood"],\n'
        '  "metadata": {"invocation": "pplacer -s raxml-info.txt -nan"},\n'
        '  "placements": [{"p": [[0, -nan, -inf], [1, nan,inf]]}]}'
    )
    parsed = parse_jplace(text)
    first, second = parsed["placements"][0]["p"]
    assert math.isnan(first[1]) and first[2] == -math.inf
    assert math.isnan(second[1]) and second[2] == math.inf
    assert parsed["metadata"]["invocation"] == "pplacer -s raxml-info.txt -nan"


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


SEPP_FIELDS = EPA_FIELDS
SEPP_TREE = (
    "((a:0.11[0],b:0.09[1])UQrYOlnDN0000020_995:0.05[2],"
    "(c:0.13[3],(d:0.08[4],e:0.15[5])UQrYOlnDN0000011_000:0.1[6]):0.0[7]);"
)


def sepp_doc(**overrides):
    return {
        "version": 1,
        "metadata": {"invocation": "run-sepp.sh"},
        "fields": SEPP_FIELDS,
        "tree": SEPP_TREE,
        "placements": [
            {
                "p": [[6, -633.4, 0.73, 0.06, 0.21], [3, -635.5, 0.09, 0.12, 0.25]],
                "nm": [["q1", 1]],
            },
            {"p": [[2, -620.0, 1.0, 0.02, 0.12]], "nm": [["q2", 1]]},
            {"p": [[7, -621.0, 1.0, 0.0, 0.12]], "nm": [["q3", 1]]},
        ],
    } | overrides


def test_from_sepp_gives_v3_with_braced_edge_numbers():
    converted = from_sepp(sepp_doc())
    assert converted["version"] == 3
    assert converted["tree"] == SEPP_TREE.replace("[", "{").replace("]", "}")
    for key in ("metadata", "fields", "placements"):
        assert converted[key] == sepp_doc()[key]


def test_from_sepp_does_not_mutate_input():
    original = sepp_doc()
    from_sepp(original)
    assert original == sepp_doc()


def test_converted_top_placements():
    result = top_placements(from_sepp(sepp_doc()))
    assert result[0] == Placement("q1", 6, 0.73, 0.06, 0.21)
    assert [p.edge_num for p in result] == [6, 2, 7]


def test_converted_doc_numbers_internal_labels():
    _, nodes = edge_numbered_tree(from_sepp(sepp_doc())["tree"])
    assert set(nodes) == set(range(8))
    assert nodes[2].label == "UQrYOlnDN0000020_995"


BACKBONE = "((a:1,b:1)X:1,(c:1,(d:1,e:1)Y:1)Z:0.5)R;"


def test_converted_doc_scores_against_truth():
    reference = "((a:1,b:1,q2:1)X:1,(c:1,(d:1,e:1,q1:1)Y:1)Z:0.5,q3:1)R;"
    rows = score_jplace(from_sepp(sepp_doc()), BACKBONE, reference, ["q1", "q2", "q3"])
    assert [r["edge_dist"] for r in rows] == [0, 0, 0]


def test_converted_root_edges_map_to_one_backbone_edge():
    mapping = jplace_edge_nodes(
        from_sepp(sepp_doc())["tree"], prepare_backbone(BACKBONE)
    )
    assert mapping[2] == mapping[7]
    assert mapping[6] != mapping[2]


@pytest.mark.parametrize(
    "bad, match",
    [
        ({"version": 3}, "version"),
        ({"version": 2}, "version"),
        ({"tree": "(a:1,b:1);"}, "numbers found"),
        ({"tree": "(a:1{0},b:1[1]);"}, "already"),
        ({"fields": ["edge_num", "like_weight_ratio"]}, "missing"),
    ],
)
def test_from_sepp_rejects(bad, match):
    with pytest.raises(ValueError, match=match):
        from_sepp(sepp_doc(**bad))
