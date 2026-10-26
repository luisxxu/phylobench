"""Stage acceptance for the full benchmark (docs/runbook-t22.md)."""

from collections import Counter

from phylobench.emp import EMP_QUERY_N, REFERENCE_ARM
from phylobench.tables import (
    Cell,
    plan_cells,
    row_count_failures,
    row_presence_failures,
    run_specs,
)

THREADS = 16
CPU_RATIO_LIMIT = 1.05
GPU_ARM = "dipper-m-gpu"
GPU_NODE = "b-gpu-21"
# apples2's guppy tog DNF is a pilot finding (docs/pilot.md), not a new failure.
ALLOWED_ERROR_DNFS = {("apples2", "convert")}


def gating_row(rows: dict, cell: Cell, arms: dict):
    return rows.get((cell.run_id, next(iter(arms[cell.arm]["phases"]))))


def coverage_failures(runs, accuracy, concordance, config: dict, arms: dict):
    """Every planned row once, and every cell scored for exactly its query_n queries."""
    failures = row_presence_failures(runs, run_specs(config, arms))
    rows = {(r["run_id"], r["phase"]): r for r in runs.to_dict("records")}
    scored = {
        "accuracy": Counter(accuracy["run_id"]),
        "concordance": Counter(concordance["run_id"]),
    }
    for cell in plan_cells(config, arms):
        gate = gating_row(rows, cell, arms)
        if gate is None:
            continue
        ok = gate["status"] == "ok"
        table = "accuracy"
        if cell.dataset in EMP_QUERY_N:
            table = "concordance"
            reference = gating_row(
                rows, cell._replace(aln="raw", arm=REFERENCE_ARM), arms
            )
            ok = ok and reference is not None and reference["status"] == "ok"
        expected = int(gate["query_n"]) if ok else 0
        failures += row_count_failures(
            cell.run_id, scored[table][cell.run_id], expected, table
        )
    return failures


def gate_failures(runs, threads: int = THREADS) -> dict[str, list[str]]:
    # evaluate runs numpy at a nominal 1 thread; its CPU time is not a tool's.
    tools = runs[runs["phase"] != "evaluate"]
    ok = runs[runs["status"] == "ok"]
    ok_tools = tools[tools["status"] == "ok"]
    place = ok[ok["phase"].str.split("+").map(lambda parts: "place" in parts)]
    gpu = place[place["arm"] == GPU_ARM]
    ratio = ok_tools["cpu_s"] / ok_tools["wall_s"]
    return {
        "mem_source": [
            f"{r.run_id} {r.phase}: mem_source {r.mem_source}"
            for r in place.itertuples()
            if r.mem_source != "cgroup"
        ],
        "cpu_ratio": [
            f"{r.run_id} {r.phase}: cpu_s/wall_s {ratio[r.Index]:.1f} > "
            f"{CPU_RATIO_LIMIT * r.threads:.1f}"
            for r in ok_tools.itertuples()
            if ratio[r.Index] > CPU_RATIO_LIMIT * r.threads
        ],
        "threads": [
            f"{r.run_id} {r.phase}: threads {r.threads}"
            for r in tools.itertuples()
            if r.threads != threads
        ],
        "gpu": [
            f"{r.run_id} {r.phase}: node {r.node}, gpu_peak_mb {r.gpu_peak_mb}"
            for r in gpu.itertuples()
            if r.node != GPU_NODE or not r.gpu_peak_mb > 0
        ],
    }


def dnf_findings(runs) -> tuple[list[str], list[str]]:
    """(error DNFs needing user acknowledgement, oom/timeout DNFs for information)."""
    errors, info = [], []
    for r in runs[runs["status"] == "DNF"].itertuples():
        line = f"{r.run_id} {r.phase}: {r.dnf_reason}"
        if r.dnf_reason != "error":
            info.append(line)
        elif (r.arm, r.phase) not in ALLOWED_ERROR_DNFS:
            errors.append(line)
    return errors, info


def stage_report(runs, accuracy, concordance, config: dict, arms: dict):
    """(lines to print, True when any hard gate failed); WARN never fails the stage."""
    results = {"coverage": coverage_failures(runs, accuracy, concordance, config, arms)}
    results.update(gate_failures(runs))
    lines = []
    for gate, failures in results.items():
        lines.append(f"{'FAIL' if failures else 'PASS'} {gate}")
        lines += [f"  {f}" for f in failures]
    errors, info = dnf_findings(runs)
    lines.append(f"{'WARN' if errors else 'PASS'} dnf_error")
    lines += [f"  {e}" for e in errors]
    lines += [f"INFO dnf {i}" for i in info]
    return lines, any(results.values())
