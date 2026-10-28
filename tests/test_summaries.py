import math

import pandas as pd
import pytest

from phylobench.summaries import (
    BIN_ORDER,
    bin_summary,
    end_to_end,
    oracle_gap,
    paired_difference,
    split_run_id,
    verdict,
)

NAN = math.nan


def acc(run_id, edge_dist, bins=None, exact=None, delta=None):
    n = len(edge_dist)
    return pd.DataFrame(
        {
            "run_id": run_id,
            "query": [f"q{i}" for i in range(n)],
            "sim": 0.9,
            "bin": bins or ["90-95"] * n,
            "edge_dist": edge_dist,
            "exact": exact if exact is not None else [float(d == 0) for d in edge_dist],
            "delta_err": delta if delta is not None else [NAN] * n,
            "label": "accuracy",
        }
    )


def cat(*frames):
    return pd.concat(frames, ignore_index=True)


def test_split_run_id_adds_five_columns():
    df = split_run_id(acc("gg2/B10k/V4-150/hmmalign/pplacer", [0.0]))
    row = df.iloc[0]
    assert (row.dataset, row.backbone, row.region, row.aln, row.arm) == (
        "gg2",
        "B10k",
        "V4-150",
        "hmmalign",
        "pplacer",
    )


def test_split_run_id_rejects_short_ids():
    with pytest.raises(ValueError):
        split_run_id(acc("gg2/B10k", [0.0]))


def test_bin_summary_never_pools_backbones():
    a = cat(
        acc("gg2/B10k/V4-150/hmmalign/x", [0.0, 0.0]),
        acc("gg2/Bfull/V4-150/hmmalign/x", [4.0, 4.0]),
    )
    out = bin_summary(a).set_index("backbone")
    assert out.loc["B10k", "mean_edge_dist"] == 0.0
    assert out.loc["Bfull", "mean_edge_dist"] == 4.0
    assert set(out["n"]) == {2}


def test_bin_summary_counts_unplaced_and_scores_them_not_exact():
    a = acc(
        "gg2/B10k/V4-150/hmmalign/x",
        [0.0, 2.0, NAN, NAN],
        exact=[1.0, 0.0, NAN, NAN],
        delta=[0.0, 2.0, NAN, NAN],
    )
    row = bin_summary(a).iloc[0]
    assert (row.n, row.n_unplaced) == (4, 2)
    assert row.mean_edge_dist == 1.0
    assert row.exact_frac == 0.25
    assert row.mean_delta_err == 1.0


def test_bin_summary_all_unplaced_has_nan_mean_and_flags_small_bins():
    row = bin_summary(acc("gg2/Bfull/V4-150/hmmalign/x", [NAN, NAN])).iloc[0]
    assert math.isnan(row.mean_edge_dist) and row.n_unplaced == 2
    assert row.small


def test_bin_summary_orders_bins_by_similarity():
    a = acc(
        "gg2/B10k/V4-150/hmmalign/x",
        [0.0] * 4,
        bins=[">=99.9", "75-80", "<75", "95-97"],
    )
    assert list(bin_summary(a)["bin"]) == ["<75", "75-80", "95-97", ">=99.9"]
    assert BIN_ORDER[0] == "<75"


def frame(values):
    return pd.DataFrame(
        {"query": [f"q{i}" for i in range(len(values))], "edge_dist": values}
    )


def test_paired_difference_constant_diff_collapses_ci():
    out = paired_difference(
        frame([3.0, 4.0, 5.0, 6.0]), frame([2.0, 3.0, 4.0, 5.0]), n_boot=200
    )
    assert out["n_pairs"] == 4
    assert out["mean"] == pytest.approx(1.0)
    assert out["ci_lo"] == pytest.approx(1.0) and out["ci_hi"] == pytest.approx(1.0)


def test_paired_difference_is_seed_deterministic():
    tool, ref = frame(list(range(0, 40, 3))), frame([1.0] * 14)
    first = paired_difference(tool, ref, n_boot=500)
    assert first == paired_difference(tool, ref, n_boot=500)
    other = paired_difference(tool, ref, n_boot=500, seed=1)
    assert (other["ci_lo"], other["ci_hi"]) != (first["ci_lo"], first["ci_hi"])


def test_paired_difference_leaves_unplaced_out_and_counts_them():
    tool = frame([1.0, NAN, 3.0, 4.0])
    ref = frame([1.0, 1.0, NAN, 2.0])
    out = paired_difference(tool, ref, n_boot=50)
    assert (out["n_pairs"], out["n_tool_unplaced"], out["n_ref_unplaced"]) == (2, 1, 1)
    assert out["mean"] == pytest.approx(1.0)


def test_paired_difference_all_unplaced_tool_has_no_pairs():
    out = paired_difference(frame([NAN, NAN]), frame([1.0, 2.0]), n_boot=50)
    assert out["n_pairs"] == 0 and out["n_tool_unplaced"] == 2
    assert math.isnan(out["mean"]) and math.isnan(out["ci_hi"])


def test_paired_difference_rejects_different_query_sets():
    tool = frame([1.0, 2.0])
    ref = frame([1.0, 2.0]).assign(query=["q0", "other"])
    with pytest.raises(ValueError, match="query"):
        paired_difference(tool, ref, n_boot=10)


def two_arm(tool_dist, ref_dist, tool="pplacer", aln="hmmalign"):
    ref = list(ref_dist)
    return cat(
        acc(f"gg2/B10k/V4-150/{aln}/{tool}", tool_dist),
        acc("gg2/B10k/V4-150/raw/sepp-prod", ref),
    )


def verdict_of(a, **kw):
    return verdict(a, **kw).set_index("arm")["verdict"].to_dict()


def test_verdict_margin_boundary():
    ref = [0.0] * 8
    assert verdict_of(two_arm([0.25] * 8, ref))["pplacer"] == "pass"
    assert verdict_of(two_arm([0.26] * 8, ref))["pplacer"] == "fail"


def test_verdict_pairs_hmmalign_challenger_with_raw_sepp_prod_and_skips_reference():
    out = verdict(two_arm([0.0] * 4, [0.0] * 4))
    assert list(out["arm"]) == ["pplacer"]
    assert out.iloc[0]["aln"] == "hmmalign"


def test_verdict_excludes_oracle():
    a = cat(
        two_arm([0.0] * 4, [0.0] * 4), acc("gg2/B10k/V4-150/oracle/pplacer", [0.0] * 4)
    )
    assert set(verdict(a)["aln"]) == {"hmmalign"}


def test_verdict_pairs_raw_challenger_on_raw():
    a = two_arm([0.0] * 4, [0.0] * 4, tool="sepp-pinned", aln="raw")
    assert verdict_of(a) == {"sepp-pinned": "pass"}


def test_verdict_unplaced_is_strict_regardless_of_ci():
    out = verdict(two_arm([0.0, 0.0, 0.0, NAN], [0.0] * 4)).iloc[0]
    assert out["verdict"] == "unplaced"
    assert out["n_tool_unplaced"] == 1 and out["n"] == 4


def test_verdict_all_unplaced_tool():
    assert verdict_of(two_arm([NAN] * 4, [0.0] * 4)) == {"pplacer": "unplaced"}


def test_verdict_no_reference_when_sepp_prod_has_unplaced():
    assert verdict_of(two_arm([0.0] * 4, [0.0, 0.0, 0.0, NAN])) == {
        "pplacer": "no-reference"
    }


def test_verdict_never_pairs_across_backbones():
    a = cat(
        acc("gg2/B10k/V4-150/hmmalign/pplacer", [0.0] * 4),
        acc("gg2/B10k/V4-150/raw/sepp-prod", [0.0] * 4),
        acc("gg2/Bfull/V4-150/hmmalign/pplacer", [9.0] * 4),
        acc("gg2/Bfull/V4-150/raw/sepp-prod", [0.0] * 4),
    )
    out = verdict(a).set_index("backbone")
    assert (
        out.loc["B10k", "verdict"] == "pass" and out.loc["Bfull", "verdict"] == "fail"
    )
    assert out.loc["Bfull", "mean"] == pytest.approx(9.0)


def test_verdict_missing_reference_cell_is_no_reference():
    assert verdict_of(acc("gg2/B10k/V4-150/hmmalign/pplacer", [0.0] * 4)) == {
        "pplacer": "no-reference"
    }


def test_oracle_gap_pairs_same_arm_oracle_minus_hmmalign():
    a = cat(
        acc("gg2/B10k/V4-150/oracle/apples2", [1.0] * 4),
        acc("gg2/B10k/V4-150/hmmalign/apples2", [3.0] * 4),
        acc("gg2/B10k/V4-150/hmmalign/pplacer", [3.0] * 4),
    )
    out = oracle_gap(a)
    assert list(out["arm"]) == ["apples2"]
    assert out.iloc[0]["mean"] == pytest.approx(-2.0)


def run_row(phase, wall, arm=None, region="V4-150", aln="hmmalign", status="ok", **kw):
    row = {
        "dataset": "gg2",
        "backbone": "B10k",
        "region": region,
        "query_n": 100,
        "tool": "t",
        "arm": arm,
        "aln": aln,
        "phase": phase,
        "wall_s": wall,
        "cpu_s": wall,
        "peak_mem_mb": 10.0,
        "mem_source": "cgroup",
        "gpu_peak_mb": NAN,
        "status": status,
        "cpu_model": "AMD",
    }
    return {**row, **kw}


def runs_of(*rows):
    return pd.DataFrame(rows)


def test_end_to_end_adds_region_align_place_and_convert():
    runs = runs_of(
        run_row("align", 5.0),
        run_row("place", 10.0, arm="pplacer"),
        run_row("convert", 1.0, arm="pplacer"),
        run_row("evaluate", 99.0, arm="pplacer"),
    )
    row = end_to_end(runs).iloc[0]
    assert row.place_s == 10.0 and row.e2e_s == 16.0
    assert row.e2e_status == "ok" and not row.convert_failed


def test_end_to_end_sepp_prod_is_its_single_phase():
    runs = runs_of(
        run_row("align", 5.0),
        run_row("align+place+convert", 50.0, arm="sepp-prod", aln="raw"),
        run_row("evaluate", 9.0, arm="sepp-prod", aln="raw"),
    )
    row = end_to_end(runs).iloc[0]
    assert row.place_s == 50.0 and row.e2e_s == 50.0


def test_end_to_end_sepp_pinned_adds_convert_but_not_hmmalign():
    runs = runs_of(
        run_row("align", 5.0),
        run_row("align+place", 40.0, arm="sepp-pinned", aln="raw"),
        run_row("convert", 2.0, arm="sepp-pinned", aln="raw"),
    )
    assert end_to_end(runs).iloc[0].e2e_s == 42.0


def test_end_to_end_dipper_has_no_convert():
    runs = runs_of(run_row("align", 5.0), run_row("place", 3.0, arm="dipper-m"))
    assert end_to_end(runs).iloc[0].e2e_s == 8.0


def test_end_to_end_failed_convert_marks_e2e_dnf_but_not_placement():
    runs = runs_of(
        run_row("align", 5.0),
        run_row("place", 10.0, arm="apples2"),
        run_row("convert", 0.4, arm="apples2", status="DNF"),
    )
    row = end_to_end(runs).iloc[0]
    assert row.e2e_status == "DNF" and row.convert_failed and row.place_status == "ok"


def test_end_to_end_failed_place_marks_dnf():
    runs = runs_of(run_row("align", 5.0), run_row("place", 10.0, arm="x", status="DNF"))
    row = end_to_end(runs).iloc[0]
    assert row.e2e_status == "DNF" and row.place_status == "DNF"


def test_end_to_end_excludes_oracle_refprep_and_hmmbuild():
    runs = runs_of(
        run_row("ref-prep", 100.0, region=None, aln=None),
        run_row("align", 100.0, region=None, aln=None),
        run_row("align", 5.0),
        run_row("place", 10.0, arm="pplacer"),
        run_row("place", 10.0, arm="pplacer", aln="oracle"),
    )
    out = end_to_end(runs)
    assert len(out) == 1 and out.iloc[0].e2e_s == 15.0


def test_end_to_end_keeps_backbones_apart():
    runs = runs_of(
        run_row("align", 5.0),
        run_row("align", 50.0, backbone="Bfull"),
        run_row("place", 10.0, arm="x"),
        run_row("place", 10.0, arm="x", backbone="Bfull"),
    )
    out = end_to_end(runs).set_index("backbone")
    assert out.loc["B10k", "e2e_s"] == 15.0 and out.loc["Bfull", "e2e_s"] == 60.0
