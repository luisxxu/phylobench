import math

import pandas as pd
import pytest
import yaml
from snakemake.exceptions import WorkflowError

from phylobench.tables import (
    ACCURACY_COLUMNS,
    CONCORDANCE_COLUMNS,
    FORCED_DNF,
    RUNS_COLUMNS,
    Cell,
    RunSpec,
    accuracy_rows,
    concordance_rows,
    forced_dnf_cell,
    load_arms,
    plan_cells,
    query_axes,
    row_count_failures,
    run_row,
    run_specs,
    score_table,
    smoke_failures,
    timing_summary,
    width_failures,
)

PPLACER = {
    "tool": "pplacer",
    "aligned": True,
    "phases": {"place": "place.run.json", "convert": "placements.tog.run.json"},
    "smoke_min_exact": 8,
}
RAW_ARM = {"tool": "rawtool", "aligned": False, "phases": {"place": "place.run.json"}}
CONFIG = {
    "datasets": ["tiny"],
    "backbones": ["B100"],
    "regions": ["full-length", "V4-150"],
    "alns": ["hmmalign", "oracle"],
    "forced_dnf": {"region": "full-length", "aln": "hmmalign"},
}
ARMS = {"pplacer": PPLACER}
QUERIES = [f"q{i}" for i in range(10)]

RUN_JSON = {
    "status": "ok",
    "wall_s": 1.5,
    "cpu_s": 1.0,
    "peak_mem_mb": 50.0,
    "mem_source": "psutil",
    "threads": 2,
    "node": "n1",
    "cpu_model": "cpu",
    "version": "1.0",
    "env_lock_sha256": "abc",
}
SPEC = RunSpec("r.json", "d/b/r/a/arm", "place", "tool", "d", "b", "r", "a", "arm")


def test_run_row_ok():
    row = run_row(SPEC, RUN_JSON, 10)
    assert row["status"] == "ok"
    assert row["dnf_reason"] == ""
    assert row["query_n"] == 10
    assert math.isnan(row["gpu_peak_mb"])


@pytest.mark.parametrize("harness_status", ["oom", "error", "timeout"])
def test_run_row_non_ok_is_dnf_with_reason(harness_status):
    row = run_row(SPEC, {**RUN_JSON, "status": harness_status}, 10)
    assert row["status"] == "DNF"
    assert row["dnf_reason"] == harness_status


def test_run_row_keeps_gpu_when_measured():
    assert run_row(SPEC, {**RUN_JSON, "gpu_peak_mb": 7.0}, 1)["gpu_peak_mb"] == 7.0


def test_run_row_columns_follow_runs_columns():
    assert tuple(run_row(SPEC, RUN_JSON, 1)) == RUNS_COLUMNS


def test_load_arms_filter_and_unknown(tmp_path):
    for name in ("b", "a"):
        (tmp_path / f"{name}.yaml").write_text(yaml.safe_dump({"tool": name}))
    assert list(load_arms(tmp_path)) == ["a", "b"]
    assert load_arms(tmp_path, only="b") == {"b": {"tool": "b"}}
    with pytest.raises(WorkflowError, match=r"unknown arm 'c'.*\['a', 'b'\]"):
        load_arms(tmp_path, only="c")


def test_plan_cells_aligned_arm_runs_over_alns_and_raw_arm_once():
    cells = plan_cells(CONFIG, {"pplacer": PPLACER, "rawtool": RAW_ARM})
    assert len(cells) == 2 * (2 + 1)
    assert Cell("tiny", "B100", "V4-150", "oracle", "pplacer") in cells
    assert Cell("tiny", "B100", "full-length", "raw", "rawtool") in cells
    assert Cell("tiny", "B100", "V4-150", "raw", "pplacer") not in cells
    assert cells[0].run_id == "tiny/B100/full-length/hmmalign/pplacer"


def test_forced_dnf_cell():
    assert forced_dnf_cell({k: v for k, v in CONFIG.items() if k != "forced_dnf"}) == []
    assert forced_dnf_cell(CONFIG) == [
        Cell("tiny", "B100", "full-length", "hmmalign", FORCED_DNF)
    ]


def test_run_specs_unique_and_complete():
    specs = run_specs(CONFIG, ARMS)
    keys = [(s.run_id, s.phase) for s in specs]
    assert len(keys) == len(set(keys))
    by_phase = {s.phase for s in specs}
    assert by_phase == {"ref-prep", "align", "place", "convert", "evaluate"}
    assert {s.tool for s in specs if s.phase == "convert"} == {"guppy"}
    assert {s.tool for s in specs if s.phase == "place"} == {"pplacer", "sleep"}
    forced = [s for s in specs if s.arm == FORCED_DNF]
    assert [s.phase for s in forced] == ["place"]
    assert forced[0].run_json.endswith("hmmalign/forced-dnf/place.run.json")


def test_accuracy_rows_join_similarity_and_keep_unplaced_nan():
    scores = [
        {"query": "a", "edge_dist": 0.0, "exact": 1.0, "delta_err": 0.0},
        {"query": "b", "edge_dist": math.nan, "exact": math.nan, "delta_err": math.nan},
    ]
    sims = pd.DataFrame(
        {"query": ["a", "b"], "sim": [0.9, 0.8], "bin": ["85-90", "80-85"]}
    )
    rows = accuracy_rows("rid", scores, sims, "accuracy")
    assert tuple(rows[0]) == ACCURACY_COLUMNS
    assert (rows[0]["sim"], rows[0]["bin"], rows[0]["label"]) == (
        0.9,
        "85-90",
        "accuracy",
    )
    assert rows[1]["sim"] == 0.8 and math.isnan(rows[1]["edge_dist"])


def test_width_failures():
    assert width_failures("f", [("a", "AC-"), ("b", "ACG")], 3) == []
    assert "b has width 2" in width_failures("f", [("a", "ACG"), ("b", "AC")], 3)[0]
    assert width_failures("f", [], 3) == ["f: empty"]
    assert width_failures("f", None, 3) == ["f: missing"]


def good_tables(config=CONFIG, arms=ARMS):
    rows = []
    for spec in run_specs(config, arms):
        forced = spec.arm == FORCED_DNF
        run = {**RUN_JSON, "status": "timeout" if forced else "ok"}
        rows.append(run_row(spec, run, 10))
    runs = pd.DataFrame(rows, columns=list(RUNS_COLUMNS))
    accuracy = pd.DataFrame(
        [
            {
                "run_id": cell.run_id,
                "query": q,
                "sim": 0.9,
                "bin": "85-90",
                "edge_dist": 0.0,
                "exact": 1.0,
                "delta_err": 0.0,
                "label": "accuracy",
            }
            for cell in plan_cells(config, arms)
            for q in QUERIES
        ],
        columns=list(ACCURACY_COLUMNS),
    )
    return runs, accuracy


def failures(runs, accuracy, arms=ARMS):
    return smoke_failures(runs, accuracy, CONFIG, arms, QUERIES)


def test_smoke_passes_on_good_tables():
    assert failures(*good_tables()) == []


def test_smoke_catches_missing_and_duplicate_and_extra_rows():
    runs, accuracy = good_tables()
    assert any("found 0" in f for f in failures(runs.iloc[1:], accuracy))
    assert any(
        "found 2" in f for f in failures(pd.concat([runs, runs.iloc[:1]]), accuracy)
    )
    extra = runs.iloc[:1].assign(phase="bogus")
    assert any("unplanned" in f for f in failures(pd.concat([runs, extra]), accuracy))


def test_smoke_catches_forced_case_not_dnf():
    runs, accuracy = good_tables()
    forced = runs["arm"] == FORCED_DNF
    assert failures(runs.assign(status=runs["status"].where(~forced, "ok")), accuracy)
    assert failures(
        runs.assign(dnf_reason=runs["dnf_reason"].where(~forced, "oom")), accuracy
    )


def test_smoke_catches_non_ok_and_null_provenance():
    runs, accuracy = good_tables()
    real = runs["arm"] == "pplacer"
    assert failures(runs.assign(status=runs["status"].where(~real, "DNF")), accuracy)
    for column in ("version", "env_lock_sha256"):
        nulled = runs.assign(**{column: runs[column].where(~real, None)})
        assert any(column in f for f in failures(nulled, accuracy))


def test_smoke_catches_accuracy_row_problems():
    runs, accuracy = good_tables()
    assert failures(runs, accuracy.iloc[1:])
    nan_dist = accuracy.assign(edge_dist=math.nan)
    assert any("NaN edge_dist" in f for f in failures(runs, nan_dist))
    forced_rows = accuracy.iloc[:1].assign(
        run_id="tiny/B100/full-length/hmmalign/forced-dnf"
    )
    assert any(
        "forced" in f for f in failures(runs, pd.concat([accuracy, forced_rows]))
    )


def test_smoke_catches_exact_below_threshold_in_full_length_only():
    runs, accuracy = good_tables()
    cell = "tiny/B100/full-length/hmmalign/pplacer"
    weak = accuracy.copy()
    weak.loc[weak["run_id"] == cell, "exact"] = [1.0] * 7 + [0.0] * 3
    assert any("exact 7 < 8" in f for f in failures(runs, weak))
    other = accuracy.copy()
    other.loc[other["run_id"] == "tiny/B100/V4-150/hmmalign/pplacer", "exact"] = 0.0
    assert failures(runs, other) == []


SEPP_PROD = {
    "tool": "sepp",
    "aligned": False,
    "phases": {"align+place+convert": "sepp.run.json"},
}
EMP_CONFIG = {
    "datasets": ["tiny", "emp-10k"],
    "backbones": ["B100"],
    "regions": ["full-length", "V4-150"],
    "alns": ["hmmalign", "oracle"],
}
EMP_ARMS = {"pplacer": PPLACER, "sepp-prod": SEPP_PROD}


def test_query_axes_emp_is_v4_hmmalign_only():
    assert query_axes(EMP_CONFIG, "emp-100k") == (["V4-150"], ["hmmalign"])
    assert query_axes(EMP_CONFIG, "tiny") == (EMP_CONFIG["regions"], EMP_CONFIG["alns"])


def test_plan_cells_emp_dataset_gets_one_region_and_aln():
    cells = plan_cells(EMP_CONFIG, EMP_ARMS)
    emp = [c for c in cells if c.dataset == "emp-10k"]
    assert [(c.region, c.aln, c.arm) for c in emp] == [
        ("V4-150", "hmmalign", "pplacer"),
        ("V4-150", "raw", "sepp-prod"),
    ]
    assert len([c for c in cells if c.dataset == "tiny"]) == 2 * 2 + 2


def test_plan_cells_non_emp_order_is_unchanged():
    assert plan_cells(CONFIG, ARMS) == [
        Cell("tiny", "B100", region, aln, "pplacer")
        for region in CONFIG["regions"]
        for aln in CONFIG["alns"]
    ]


def test_emp_run_specs_skip_refprep_and_hmmbuild_and_use_concordance():
    specs = [s for s in run_specs(EMP_CONFIG, EMP_ARMS) if s.dataset == "emp-10k"]
    assert {s.phase for s in specs} == {
        "align",
        "place",
        "convert",
        "align+place+convert",
        "evaluate",
    }
    assert [s.run_json for s in specs if s.phase == "align"] == [
        "results/emp-10k/B100/V4-150/align/hmmalign.run.json"
    ]
    evaluate = {s.run_json: s for s in specs if s.phase == "evaluate"}
    assert set(evaluate) == {
        "results/emp-10k/B100/V4-150/hmmalign/pplacer/concordance.run.json",
        "results/emp-10k/B100/V4-150/raw/sepp-prod/concordance.run.json",
    }
    assert {s.tool for s in evaluate.values()} == {"phylobench"}


def test_non_emp_specs_still_use_accuracy_evaluate():
    specs = run_specs(EMP_CONFIG, EMP_ARMS)
    tiny = [s for s in specs if s.dataset == "tiny" and s.phase == "evaluate"]
    assert all(s.run_json.endswith("/evaluate.run.json") for s in tiny)
    assert any(s.phase == "ref-prep" and s.dataset == "tiny" for s in specs)


def test_emp_without_sepp_prod_raises():
    with pytest.raises(WorkflowError, match="sepp-prod"):
        run_specs(EMP_CONFIG, ARMS)


def test_emp_single_arm_run_may_omit_sepp_prod():
    assert run_specs({**EMP_CONFIG, "arm": "pplacer"}, ARMS)


def test_gg2_without_sepp_prod_is_fine():
    assert run_specs(CONFIG, ARMS)


def test_score_table_follows_the_evaluate_run_json():
    emp = RunSpec("d/c/concordance.run.json", "r", "evaluate", "t", *"dbraa")
    acc = RunSpec("d/c/evaluate.run.json", "r", "evaluate", "t", *"dbraa")
    assert score_table(emp) == "d/c/concordance.tsv"
    assert score_table(acc) == "d/c/accuracy.tsv"


def test_concordance_rows_label_and_nan():
    scores = [
        {"query": "a", "edge_dist": 2.0},
        {"query": "b", "edge_dist": math.nan},
    ]
    sims = pd.DataFrame({"query": ["a", "b"], "sim": [0.9, 0.8], "bin": ["a", "b"]})
    rows = concordance_rows("rid", "ref", scores, sims)
    assert tuple(rows[0]) == CONCORDANCE_COLUMNS
    assert {r["label"] for r in rows} == {"concordance"}
    assert rows[0]["reference_run_id"] == "ref" and rows[0]["sim"] == 0.9
    assert math.isnan(rows[1]["edge_dist"])


def timing_runs(*rows):
    base = {
        "dataset": "gg2",
        "backbone": "B10k",
        "region": "V4-150",
        "aln": "hmmalign",
        "arm": "pplacer",
        "phase": "place",
        "status": "ok",
        "peak_mem_mb": 100.0,
    }
    return pd.DataFrame([base | row for row in rows])


def test_timing_summary_mean_and_cv_over_ok_rows():
    runs = timing_runs({"wall_s": 10.0}, {"wall_s": 12.0}, {"wall_s": 14.0})
    row = timing_summary(runs).iloc[0]
    assert (row["n"], row["n_ok"], row["wall_mean"]) == (3, 3, 12.0)
    assert row["wall_cv"] == pytest.approx(1 / 6)
    assert row["peak_mem_max"] == 100.0


def test_timing_summary_counts_dnf_in_n_not_n_ok():
    runs = timing_runs(
        {"wall_s": 10.0},
        {"wall_s": 14.0},
        {"wall_s": 99.0, "status": "DNF", "peak_mem_mb": 300.0},
    )
    row = timing_summary(runs).iloc[0]
    assert (row["n"], row["n_ok"], row["wall_mean"]) == (3, 2, 12.0)
    assert row["peak_mem_max"] == 300.0


def test_timing_summary_cv_is_nan_below_two_ok_rows():
    runs = timing_runs({"wall_s": 10.0}, {"wall_s": 5.0, "status": "DNF"})
    row = timing_summary(runs).iloc[0]
    assert row["n_ok"] == 1 and math.isnan(row["wall_cv"])


def test_timing_summary_keeps_all_dnf_cell_and_drops_other_phases():
    runs = timing_runs(
        {"wall_s": 5.0, "status": "DNF"},
        {"wall_s": 1.0, "phase": "evaluate", "arm": "epa-ng"},
    )
    summary = timing_summary(runs)
    assert list(summary["arm"]) == ["pplacer"]
    assert summary.iloc[0]["n_ok"] == 0 and math.isnan(summary.iloc[0]["wall_mean"])


def test_timing_summary_includes_combined_place_phases():
    runs = timing_runs(
        {"wall_s": 5.0, "phase": "align+place+convert", "arm": "sepp-prod"},
        {"wall_s": 6.0, "phase": "align+place", "arm": "sepp-pinned"},
        {"wall_s": 1.0, "phase": "convert", "arm": "sepp-pinned"},
    )
    summary = timing_summary(runs).set_index("arm")
    assert summary.loc["sepp-prod", "wall_mean"] == 5.0
    assert summary.loc["sepp-pinned", "wall_mean"] == 6.0


TWO_ARMS = {"pplacer": PPLACER, "rawtool": RAW_ARM}
TWO_BACKBONES = {
    **CONFIG,
    "backbones": ["B10k", "Bfull"],
    "regions": ["full-length", "V4-150"],
}


def planned_with(exclude, config=TWO_BACKBONES):
    return plan_cells({**config, "exclude": exclude}, TWO_ARMS)


def test_exclude_ands_keys_within_an_entry():
    kept = planned_with([{"backbone": ["Bfull"], "arm": ["rawtool"]}])
    assert Cell("tiny", "B10k", "V4-150", "raw", "rawtool") in kept
    assert Cell("tiny", "Bfull", "V4-150", "raw", "rawtool") not in kept
    assert Cell("tiny", "Bfull", "V4-150", "oracle", "pplacer") in kept


def test_exclude_ors_values_within_a_key_and_across_entries():
    kept = planned_with(
        [{"region": ["full-length", "V4-150"], "arm": ["rawtool"]}, {"aln": ["oracle"]}]
    )
    assert {(c.arm, c.aln) for c in kept} == {("pplacer", "hmmalign")}


def test_exclude_absent_changes_nothing():
    assert planned_with([]) == plan_cells(TWO_BACKBONES, TWO_ARMS)


def test_excluded_backbone_gets_no_shared_specs():
    config = {**TWO_BACKBONES, "exclude": [{"backbone": ["Bfull"]}]}
    specs = run_specs(config, TWO_ARMS)
    assert {s.backbone for s in specs} == {"B10k"}


def test_hmmalign_spec_only_for_planned_regions_and_kept_for_raw_only_arms():
    config = {**TWO_BACKBONES, "exclude": [{"region": ["V4-150"], "arm": ["pplacer"]}]}
    aligned = [s for s in run_specs(config, TWO_ARMS) if s.run_id.endswith("hmmalign")]
    assert {s.region for s in aligned} == {"full-length", "V4-150"}
    only_raw = {**TWO_BACKBONES, "exclude": [{"arm": ["pplacer"]}]}
    assert {
        s.region for s in run_specs(only_raw, TWO_ARMS) if s.run_id.endswith("hmmalign")
    } == {"full-length", "V4-150"}


def test_run_specs_without_exclude_keep_shared_specs_first_in_config_order():
    specs = run_specs(TWO_BACKBONES, TWO_ARMS)
    assert [(s.backbone, s.phase) for s in specs[:3]] == [("B10k", "ref-prep")] + [
        ("B10k", "align")
    ] * 2


def test_row_count_failures():
    assert row_count_failures("r", 3, 3, "accuracy") == []
    assert row_count_failures("r", 2, 3, "accuracy") == [
        "r: 2 accuracy rows, expected 3"
    ]
