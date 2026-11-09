import math

import pandas as pd
import pytest

from phylobench.summaries import within_k_difference

NAN = math.nan
REF_CELL = ("gg2-wol", "B10k", "full-length", "raw", "sepp-prod")
TOOL_CELL = ("gg2-wol", "B10k", "full-length", "hmmalign", "bscampp")


def cell(arm, edge_dist, *, dataset="gg2-wol", backbone="B10k"):
    aln = "raw" if arm == "sepp-prod" else "hmmalign"
    return pd.DataFrame(
        {
            "run_id": f"{dataset}/{backbone}/full-length/{aln}/{arm}",
            "query": [f"q{i}" for i in range(len(edge_dist))],
            "edge_dist": edge_dist,
        }
    )


def pair(ref, tool):
    return pd.concat([cell("sepp-prod", ref), cell("bscampp", tool)], ignore_index=True)


def groups_of(n, size=1):
    return pd.Series({f"q{i}": f"g{i // size}" for i in range(n)})


def test_identical_arms_give_zero_difference_and_bound():
    ed = [0.0, 1.0, 5.0, 2.0] * 5
    row = within_k_difference(pair(ed, ed), groups_of(20)).iloc[0]
    assert (row["diff"], row["lower"]) == (0, 0)
    assert (row["n"], row["n_clusters"]) == (20, 20)
    assert row["verdict"] == "non-inferior"


def test_known_difference_and_verdict():
    tool = [0.0] * 8 + [9.0, 9.0]
    row = within_k_difference(pair([0.0] * 10, tool), groups_of(10)).iloc[0]
    assert row["diff"] == pytest.approx(-0.2)
    assert row["lower"] < -0.02
    assert row["verdict"] == "not-shown"


def test_edge_dist_equal_to_k_is_a_hit():
    row = within_k_difference(pair([3.0] * 4, [4.0] * 4), groups_of(4), k=3).iloc[0]
    assert row["diff"] == -1.0


def test_unplaced_challenger_query_is_a_miss():
    row = within_k_difference(pair([0.0] * 4, [0.0, 0.0, 0.0, NAN]), groups_of(4)).iloc[
        0
    ]
    assert row["diff"] == -0.25


def test_planned_challenger_without_rows_is_all_miss():
    row = within_k_difference(
        cell("sepp-prod", [0.0] * 4),
        groups_of(4),
        planned_cells=[REF_CELL, TOOL_CELL],
    ).iloc[0]
    assert (row["arm"], row["diff"], row["n"]) == ("bscampp", -1.0, 4)


def test_unplaced_reference_query_is_a_miss():
    row = within_k_difference(pair([0.0, NAN], [0.0, 0.0]), groups_of(2)).iloc[0]
    assert row["diff"] == 0.5
    assert row["verdict"] != "undecidable"


def test_dnf_reference_is_undecidable():
    out = within_k_difference(
        cell("bscampp", [0.0, 0.0]),
        groups_of(2),
        planned_cells=[REF_CELL, TOOL_CELL],
    )
    assert out["verdict"].tolist() == ["undecidable"]


def test_clustering_widens_the_interval():
    tool = [9.0 if i // 10 < 5 else 0.0 for i in range(200)]
    accuracy = pair([0.0] * 200, tool)
    clustered = within_k_difference(accuracy, groups_of(200, 10)).iloc[0]
    per_query = within_k_difference(accuracy, groups_of(200)).iloc[0]
    assert clustered["n_clusters"] == 20
    assert clustered["diff"] == per_query["diff"] == -0.25
    assert clustered["lower"] < per_query["lower"]


def test_seed_determines_the_bound():
    tool = [0.0, 9.0, 1.0, 9.0, 0.0, 2.0] * 50
    accuracy = pair([0.0] * 300, tool)
    args = (accuracy, groups_of(300))
    a = within_k_difference(*args, n_boot=500, seed=1)
    assert a.equals(within_k_difference(*args, n_boot=500, seed=1))
    other = within_k_difference(*args, n_boot=500, seed=2)
    assert a["lower"].iloc[0] != other["lower"].iloc[0]


def test_never_pools_backbones_or_datasets():
    frames = []
    for dataset, backbone in [
        ("gg2-wol", "B10k"),
        ("gg2-wol", "Bfull"),
        ("gg2", "B10k"),
    ]:
        for arm in ("sepp-prod", "bscampp"):
            frames.append(cell(arm, [0.0] * 4, dataset=dataset, backbone=backbone))
    out = within_k_difference(pd.concat(frames, ignore_index=True), groups_of(4))
    assert len(out) == 3
    assert set(out["n"]) == {4}


def test_query_without_a_group_is_an_error():
    with pytest.raises(ValueError, match="group"):
        within_k_difference(pair([0.0] * 3, [0.0] * 3), groups_of(2))
