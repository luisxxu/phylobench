import math

import pandas as pd
import pytest
import yaml
from snakemake.exceptions import WorkflowError

from phylobench.tables import (
    ACCURACY_COLUMNS,
    FORCED_DNF,
    RUNS_COLUMNS,
    Cell,
    RunSpec,
    accuracy_rows,
    forced_dnf_cell,
    load_arms,
    plan_cells,
    run_row,
    run_specs,
    smoke_failures,
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
