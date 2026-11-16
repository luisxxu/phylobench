from pathlib import Path

import pytest

from phylobench.regions import read_fasta, upp_query_rows

FIXTURE = Path(__file__).parents[1] / "tests" / "fixtures" / "tiny"
BACKBONE = {"b1": "ACGT-a", "b2": "AC-TTG"}
MASKED = {"b1": "ACGT-A", "b2": "AC-TTG", "q1": "A-GTTG", "q2": "ACG--T"}


def test_returns_present_queries_in_query_order():
    rows = upp_query_rows(MASKED, BACKBONE, ["q2", "q1"])
    assert rows == [("q2", "ACG--T"), ("q1", "A-GTTG")]


def test_lowercase_backbone_matches_its_uppercased_row():
    assert upp_query_rows({"b1": "ACGT-A", "b2": "AC-TTG"}, BACKBONE, []) == []


def test_missing_query_is_omitted():
    masked = {k: v for k, v in MASKED.items() if k != "q1"}
    assert upp_query_rows(masked, BACKBONE, ["q1", "q2"]) == [("q2", "ACG--T")]


def test_width_mismatch_raises():
    with pytest.raises(ValueError, match="q1.*width"):
        upp_query_rows({**MASKED, "q1": "A-GTT"}, BACKBONE, ["q1", "q2"])


def test_changed_backbone_row_raises():
    with pytest.raises(ValueError, match="b2"):
        upp_query_rows({**MASKED, "b2": "AC-TTA"}, BACKBONE, ["q1", "q2"])


def test_shifted_backbone_row_raises():
    with pytest.raises(ValueError, match="b2"):
        upp_query_rows({**MASKED, "b2": "-AC-TT"}, BACKBONE, ["q1", "q2"])


def test_missing_backbone_row_raises():
    masked = {k: v for k, v in MASKED.items() if k != "b1"}
    with pytest.raises(ValueError, match="b1"):
        upp_query_rows(masked, BACKBONE, ["q1", "q2"])


@pytest.mark.parametrize("bad", ["A-gTTG", "A.GTTG"])
def test_lowercase_or_dot_in_a_query_row_raises(bad):
    with pytest.raises(ValueError, match="q1"):
        upp_query_rows({**MASKED, "q1": bad}, BACKBONE, ["q1", "q2"])


def test_unknown_name_raises():
    with pytest.raises(ValueError, match="stray"):
        upp_query_rows({**MASKED, "stray": "ACGT-A"}, BACKBONE, ["q1", "q2"])


def test_tiny_fixture_backbone_round_trips():
    backbone = dict(read_fasta(FIXTURE / "aligned-dna-sequences.fasta"))
    queries = [n for n, _ in read_fasta(FIXTURE / "queries.fasta")]
    masked = {n: r.upper() for n, r in backbone.items()}
    assert upp_query_rows(masked, backbone, queries) == []
