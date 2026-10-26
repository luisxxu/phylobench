import pandas as pd
import pytest

from phylobench.gates import (
    coverage_failures,
    dnf_findings,
    gate_failures,
    stage_report,
)
from phylobench.tables import (
    ACCURACY_COLUMNS,
    CONCORDANCE_COLUMNS,
    RUNS_COLUMNS,
    run_row,
    run_specs,
)

PHASES = {"place": "place.run.json", "convert": "placements.tog.run.json"}
ARMS = {
    "pplacer": {"tool": "pplacer", "aligned": True, "phases": PHASES},
    "apples2": {"tool": "apples2", "aligned": True, "phases": PHASES},
    "dipper-m-gpu": {
        "tool": "dipper",
        "aligned": True,
        "phases": {"place": "place.run.json"},
    },
    "sepp-prod": {
        "tool": "sepp",
        "aligned": False,
        "phases": {"align+place+convert": "sepp.run.json"},
    },
}
CONFIG = {
    "datasets": ["gg2", "emp-10k"],
    "backbones": ["B10k"],
    "regions": ["V4-150"],
    "alns": ["hmmalign"],
}
QUERY_N = 4
RUN_JSON = {
    "status": "ok",
    "wall_s": 10.0,
    "cpu_s": 100.0,
    "peak_mem_mb": 5.0,
    "mem_source": "cgroup",
    "threads": 16,
    "node": "b-cn-1",
    "cpu_model": "cpu",
    "version": "v",
    "env_lock_sha256": "e",
}


def good_runs(overrides=None):
    overrides = overrides or {}
    rows = []
    for spec in run_specs(CONFIG, ARMS):
        run = {**RUN_JSON}
        if spec.arm == "dipper-m-gpu":
            run |= {"node": "b-gpu-21", "gpu_peak_mb": 900.0}
        run |= overrides.get((spec.run_id, spec.phase), {})
        if spec.phase == "evaluate":
            run["threads"] = 1
        rows.append(run_row(spec, run, QUERY_N))
    return pd.DataFrame(rows, columns=list(RUNS_COLUMNS))


def score_tables(runs, drop=()):
    def table(columns, run_ids):
        rows = [
            {"run_id": r, "query": f"q{i}"} for r in run_ids for i in range(QUERY_N)
        ]
        return pd.DataFrame(rows, columns=list(columns))

    def scored(dataset_is_emp):
        return [
            run_id
            for run_id in runs.loc[runs["phase"] == "evaluate", "run_id"]
            if run_id not in drop and ("/emp-" in f"/{run_id}") == dataset_is_emp
        ]

    return table(ACCURACY_COLUMNS, scored(False)), table(
        CONCORDANCE_COLUMNS, scored(True)
    )


def coverage(runs, accuracy, concordance):
    return coverage_failures(runs, accuracy, concordance, CONFIG, ARMS)


def test_coverage_passes_on_complete_tables():
    runs = good_runs()
    assert coverage(runs, *score_tables(runs)) == []


def test_coverage_flags_missing_duplicate_and_unplanned_rows():
    runs = good_runs()
    accuracy, concordance = score_tables(runs)
    assert any("found 0" in f for f in coverage(runs.iloc[1:], accuracy, concordance))
    assert any(
        "found 2" in f
        for f in coverage(pd.concat([runs, runs.iloc[:1]]), accuracy, concordance)
    )
    extra = runs.iloc[:1].assign(phase="bogus")
    assert any(
        "unplanned" in f
        for f in coverage(pd.concat([runs, extra]), accuracy, concordance)
    )


def test_coverage_flags_short_score_tables():
    runs = good_runs()
    accuracy, concordance = score_tables(runs)
    assert any(
        "accuracy rows, expected 4" in f
        for f in coverage(runs, accuracy.iloc[1:], concordance)
    )
    assert any(
        "concordance rows, expected 4" in f
        for f in coverage(runs, accuracy, concordance.iloc[1:])
    )


def test_coverage_dnf_gate_must_have_no_score_rows():
    gate = "gg2/B10k/V4-150/hmmalign/pplacer", "place"
    runs = good_runs({gate: {"status": "oom"}})
    accuracy, concordance = score_tables(runs)
    assert any(
        "pplacer: 4 accuracy rows, expected 0" in f
        for f in coverage(runs, accuracy, concordance)
    )
    accuracy, concordance = score_tables(runs, drop={gate[0]})
    assert coverage(runs, accuracy, concordance) == []


def test_coverage_emp_cell_needs_sepp_prod_ok_too():
    reference = "emp-10k/B10k/V4-150/raw/sepp-prod", "align+place+convert"
    runs = good_runs({reference: {"status": "timeout"}})
    accuracy, concordance = score_tables(
        runs, drop={"emp-10k/B10k/V4-150/hmmalign/pplacer"}
    )
    problems = coverage(runs, accuracy, concordance)
    assert [p for p in problems if "pplacer" in p] == []
    assert any("sepp-prod: 4 concordance rows, expected 0" in p for p in problems)


def test_gates_pass_on_good_runs():
    assert not any(gate_failures(good_runs()).values())


@pytest.mark.parametrize(
    ("gate", "override"),
    [
        ("mem_source", {"mem_source": "psutil"}),
        ("cpu_ratio", {"cpu_s": 17.0 * 10}),
        ("threads", {"threads": 8, "cpu_s": 10.0}),
        ("gpu", {"node": "b-gpu-22"}),
        ("gpu", {"gpu_peak_mb": 0.0}),
    ],
)
def test_each_gate_fails_on_its_violation(gate, override):
    key = {
        "mem_source": ("gg2/B10k/V4-150/hmmalign/pplacer", "place"),
        "cpu_ratio": ("gg2/B10k/V4-150/hmmalign/pplacer", "place"),
        "threads": ("gg2/B10k/V4-150/hmmalign/pplacer", "convert"),
        "gpu": ("gg2/B10k/V4-150/hmmalign/dipper-m-gpu", "place"),
    }[gate]
    failures = gate_failures(good_runs({key: override}))
    assert len(failures[gate]) == 1 and key[0] in failures[gate][0]
    assert sum(map(bool, failures.values())) == 1


def test_cpu_ratio_boundary_and_evaluate_threads_exempt():
    key = ("gg2/B10k/V4-150/hmmalign/pplacer", "place")
    assert not gate_failures(good_runs({key: {"cpu_s": 168.0}}))["cpu_ratio"]
    assert not gate_failures(good_runs())["threads"]


def test_mem_source_only_gates_ok_place_rows():
    key = ("gg2/B10k/V4-150/hmmalign/pplacer", "place")
    runs = good_runs({key: {"mem_source": "psutil", "status": "oom"}})
    assert not gate_failures(runs)["mem_source"]


def test_dnf_findings_split_error_from_oom_timeout_and_allow_apples2_convert():
    runs = good_runs(
        {
            ("gg2/B10k/V4-150/hmmalign/pplacer", "place"): {"status": "error"},
            ("gg2/B10k/V4-150/hmmalign/apples2", "convert"): {"status": "error"},
            ("gg2/B10k/V4-150/hmmalign/apples2", "place"): {"status": "oom"},
            ("emp-10k/B10k/V4-150/hmmalign/pplacer", "place"): {"status": "timeout"},
        }
    )
    errors, info = dnf_findings(runs)
    assert errors == ["gg2/B10k/V4-150/hmmalign/pplacer place: error"]
    assert sorted(info) == [
        "emp-10k/B10k/V4-150/hmmalign/pplacer place: timeout",
        "gg2/B10k/V4-150/hmmalign/apples2 place: oom",
    ]


def test_stage_report_warn_does_not_fail_but_fail_does():
    runs = good_runs(
        {("gg2/B10k/V4-150/hmmalign/pplacer", "place"): {"status": "error"}}
    )
    accuracy, concordance = score_tables(
        runs, drop={"gg2/B10k/V4-150/hmmalign/pplacer"}
    )
    lines, failed = stage_report(runs, accuracy, concordance, CONFIG, ARMS)
    assert not failed and "WARN dnf_error" in lines and "PASS coverage" in lines
    lines, failed = stage_report(runs, accuracy.iloc[1:], concordance, CONFIG, ARMS)
    assert failed and "FAIL coverage" in lines
