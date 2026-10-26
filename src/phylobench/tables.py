"""Run planning and the tidy runs/accuracy tables (docs/harness.md section 7)."""

import math
from collections import Counter
from itertools import product
from pathlib import Path
from typing import NamedTuple

import pandas as pd
import yaml
from snakemake.exceptions import WorkflowError

from phylobench.emp import EMP_QUERY_N, EMP_REGION, REFERENCE_ARM

RUNS_COLUMNS = (
    "dataset",
    "backbone",
    "region",
    "query_n",
    "tool",
    "arm",
    "aln",
    "phase",
    "wall_s",
    "cpu_s",
    "peak_mem_mb",
    "mem_source",
    "gpu_peak_mb",
    "status",
    "threads",
    "node",
    "cpu_model",
    "version",
    "run_id",
    "dnf_reason",
    "env_lock_sha256",
)
ACCURACY_COLUMNS = (
    "run_id",
    "query",
    "sim",
    "bin",
    "edge_dist",
    "exact",
    "delta_err",
    "label",
)
CONCORDANCE_COLUMNS = (
    "run_id",
    "reference_run_id",
    "query",
    "sim",
    "bin",
    "edge_dist",
    "label",
)
CONCORDANCE_RUN = "concordance.run.json"
FORCED_DNF = "forced-dnf"
RAW = "raw"
FULL_LENGTH = "full-length"
TIMING_KEYS = ["dataset", "backbone", "region", "aln", "arm"]


class Cell(NamedTuple):
    dataset: str
    backbone: str
    region: str
    aln: str
    arm: str

    @property
    def run_id(self) -> str:
        return "/".join(self)

    @property
    def dir(self) -> str:
        return f"results/{self.run_id}"


class RunSpec(NamedTuple):
    run_json: str
    run_id: str
    phase: str
    tool: str
    dataset: str
    backbone: str
    region: str
    aln: str
    arm: str


def load_arms(arms_dir, only: str | None = None) -> dict[str, dict]:
    arms = {
        p.stem: yaml.safe_load(p.read_text())
        for p in sorted(Path(arms_dir).glob("*.yaml"))
    }
    if only is not None and only not in arms:
        raise WorkflowError(f"unknown arm {only!r}; available: {sorted(arms)}")
    return arms if only is None else {only: arms[only]}


def query_axes(config: dict, dataset: str) -> tuple[list[str], list[str]]:
    """(regions, alns) a dataset runs over; EMP has one region and no oracle."""
    if dataset in EMP_QUERY_N:
        return [EMP_REGION], ["hmmalign"]
    return config["regions"], config["alns"]


def excluded(cell: Cell, patterns: list[dict]) -> bool:
    """True when every key of some pattern lists the cell's value."""
    return any(
        all(getattr(cell, key) in values for key, values in pattern.items())
        for pattern in patterns
    )


def plan_cells(config: dict, arms: dict[str, dict]) -> list[Cell]:
    patterns = config.get("exclude", [])
    cells = []
    for dataset, backbone in product(config["datasets"], config["backbones"]):
        regions, alns = query_axes(config, dataset)
        for region, (name, arm) in product(regions, arms.items()):
            for aln in alns if arm["aligned"] else [RAW]:
                cell = Cell(dataset, backbone, region, aln, name)
                if not excluded(cell, patterns):
                    cells.append(cell)
    return cells


def forced_dnf_cell(config: dict) -> list[Cell]:
    forced = config.get("forced_dnf")
    if not forced:
        return []
    return [
        Cell(
            config["datasets"][0],
            config["backbones"][0],
            forced["region"],
            forced["aln"],
            FORCED_DNF,
        )
    ]


def run_specs(config: dict, arms: dict[str, dict]) -> list[RunSpec]:
    # --config arm=X plans one arm; the concordance rule still builds sepp-prod itself.
    if (
        any(d in EMP_QUERY_N for d in config["datasets"])
        and REFERENCE_ARM not in arms
        and "arm" not in config
    ):
        raise WorkflowError(
            f"EMP datasets need the {REFERENCE_ARM!r} arm as the concordance reference"
        )
    specs = []
    cells = plan_cells(config, arms)
    for dataset, backbone in product(config["datasets"], config["backbones"]):
        base = f"results/{dataset}/{backbone}"
        shared = f"{dataset}/{backbone}"
        planned = {
            c.region for c in cells if (c.dataset, c.backbone) == (dataset, backbone)
        }
        if not planned:
            continue
        regions = [r for r in query_axes(config, dataset)[0] if r in planned]
        if dataset not in EMP_QUERY_N:
            specs.append(
                RunSpec(
                    f"results/ref-prep/{shared}/me-tree.run.json",
                    shared,
                    "ref-prep",
                    "fasttree",
                    dataset,
                    backbone,
                    "",
                    "",
                    "",
                )
            )
            specs.append(
                RunSpec(
                    f"{base}/align/hmmbuild.run.json",
                    shared,
                    "align",
                    "hmmer",
                    dataset,
                    backbone,
                    "",
                    "",
                    "",
                )
            )
        for region in regions:
            specs.append(
                RunSpec(
                    f"{base}/{region}/align/hmmalign.run.json",
                    f"{shared}/{region}/hmmalign",
                    "align",
                    "hmmer",
                    dataset,
                    backbone,
                    region,
                    "hmmalign",
                    "",
                )
            )
    for cell in cells:
        arm = arms[cell.arm]
        for phase, filename in arm["phases"].items():
            tool = "guppy" if phase == "convert" else arm["tool"]
            specs.append(cell_spec(cell, f"{cell.dir}/{filename}", phase, tool))
        scorer = "concordance" if cell.dataset in EMP_QUERY_N else "evaluate"
        specs.append(
            cell_spec(cell, f"{cell.dir}/{scorer}.run.json", "evaluate", "phylobench")
        )
    for cell in forced_dnf_cell(config):
        specs.append(cell_spec(cell, f"{cell.dir}/place.run.json", "place", "sleep"))
    return specs


def score_table(spec: RunSpec) -> str:
    """The per-cell table written next to an evaluate-phase run.json."""
    name = (
        "concordance.tsv"
        if Path(spec.run_json).name == CONCORDANCE_RUN
        else "accuracy.tsv"
    )
    return str(Path(spec.run_json).with_name(name))


def cell_spec(cell: Cell, run_json: str, phase: str, tool: str) -> RunSpec:
    return RunSpec(run_json, cell.run_id, phase, tool, *cell)


def run_row(spec: RunSpec, run: dict, query_n) -> dict:
    ok = run["status"] == "ok"
    return {
        "dataset": spec.dataset,
        "backbone": spec.backbone,
        "region": spec.region,
        "query_n": query_n,
        "tool": spec.tool,
        "arm": spec.arm,
        "aln": spec.aln,
        "phase": spec.phase,
        "wall_s": run["wall_s"],
        "cpu_s": run["cpu_s"],
        "peak_mem_mb": run["peak_mem_mb"],
        "mem_source": run["mem_source"],
        "gpu_peak_mb": run.get("gpu_peak_mb", math.nan),
        "status": "ok" if ok else "DNF",
        "threads": run["threads"],
        "node": run["node"],
        "cpu_model": run["cpu_model"],
        "version": run["version"],
        "run_id": spec.run_id,
        "dnf_reason": "" if ok else run["status"],
        "env_lock_sha256": run["env_lock_sha256"],
    }


def accuracy_rows(run_id: str, scores: list[dict], sims, label: str) -> list[dict]:
    """Join per-query scores with the similarity table (query, sim, bin)."""
    sim_of = dict(zip(sims["query"], sims["sim"], strict=True))
    bin_of = dict(zip(sims["query"], sims["bin"], strict=True))
    return [
        {
            "run_id": run_id,
            "query": s["query"],
            "sim": sim_of.get(s["query"], math.nan),
            "bin": bin_of.get(s["query"], ""),
            "edge_dist": s["edge_dist"],
            "exact": s["exact"],
            "delta_err": s["delta_err"],
            "label": label,
        }
        for s in scores
    ]


def concordance_rows(
    run_id: str, reference_run_id: str, scores: list[dict], sims
) -> list[dict]:
    """Join per-query edge distances to the reference placement with similarity."""
    sim_of = dict(zip(sims["query"], sims["sim"], strict=True))
    bin_of = dict(zip(sims["query"], sims["bin"], strict=True))
    return [
        {
            "run_id": run_id,
            "reference_run_id": reference_run_id,
            "query": s["query"],
            "sim": sim_of.get(s["query"], math.nan),
            "bin": bin_of.get(s["query"], ""),
            "edge_dist": s["edge_dist"],
            "label": "concordance",
        }
        for s in scores
    ]


def width_failures(name: str, records, width: int) -> list[str]:
    """records: (sequence name, row) pairs, or None when the file is missing."""
    if records is None:
        return [f"{name}: missing"]
    if not records:
        return [f"{name}: empty"]
    return [
        f"{name}: {seq_name} has width {len(row)}, backbone has {width}"
        for seq_name, row in records
        if len(row) != width
    ]


def row_count_failures(run_id: str, found: int, expected: int, table: str) -> list[str]:
    if found == expected:
        return []
    return [f"{run_id}: {found} {table} rows, expected {expected}"]


def row_presence_failures(runs, specs: list[RunSpec]) -> list[str]:
    """Every planned (run_id, phase) exactly once in runs, and nothing else."""
    failures = []
    planned = {(s.run_id, s.phase) for s in specs}
    seen = Counter(zip(runs["run_id"], runs["phase"], strict=True))
    for key in sorted(planned):
        if seen[key] != 1:
            failures.append(f"{key}: expected exactly one row, found {seen[key]}")
    for key in sorted(set(seen) - planned):
        failures.append(f"{key}: unplanned row")
    return failures


def smoke_failures(runs, accuracy, config: dict, arms: dict, queries: list[str]):
    failures = row_presence_failures(runs, run_specs(config, arms))

    forced_ids = {c.run_id for c in forced_dnf_cell(config)}
    for row in runs.to_dict("records"):
        key = f"{row['run_id']} {row['phase']}"
        if row["run_id"] in forced_ids:
            if row["status"] != "DNF" or row["dnf_reason"] != "timeout":
                failures.append(
                    f"{key}: forced case must be DNF/timeout, got "
                    f"{row['status']}/{row['dnf_reason']}"
                )
            continue
        if row["status"] != "ok":
            failures.append(f"{key}: status {row['status']} ({row['dnf_reason']})")
        for column in ("version", "env_lock_sha256"):
            if not isinstance(row[column], str) or not row[column]:
                failures.append(f"{key}: {column} is null")

    wanted = sorted(queries)
    for cell in plan_cells(config, arms):
        scored = accuracy[accuracy["run_id"] == cell.run_id]
        count = row_count_failures(cell.run_id, len(scored), len(wanted), "accuracy")
        if count:
            failures += count
        elif sorted(scored["query"]) != wanted:
            failures.append(f"{cell.run_id}: queries differ from the fixture queries")
        elif scored["edge_dist"].isna().any():
            failures.append(f"{cell.run_id}: NaN edge_dist")
    for run_id in sorted(forced_ids & set(accuracy["run_id"])):
        failures.append(f"{run_id}: forced case must have no accuracy rows")

    for cell in plan_cells(config, arms):
        floor = arms[cell.arm].get("smoke_min_exact")
        if floor is None or cell.region != FULL_LENGTH:
            continue
        exact = accuracy.loc[accuracy["run_id"] == cell.run_id, "exact"].sum()
        if exact < floor:
            failures.append(f"{cell.run_id}: exact {exact:g} < {floor}")
    return failures


def timing_summary(runs: pd.DataFrame) -> pd.DataFrame:
    """Place-phase wall time per cell over repeated runs; CV is over ok rows only.

    SEPP arms measure placement inside one combined phase (e.g. `align+place+convert`).
    """
    place = runs[runs["phase"].str.split("+").map(lambda parts: "place" in parts)]
    grouped = place.groupby(TIMING_KEYS)
    ok = place[place["status"] == "ok"].groupby(TIMING_KEYS)["wall_s"]
    summary = pd.DataFrame(
        {
            "n": grouped.size(),
            "n_ok": ok.size(),
            "wall_mean": ok.mean(),
            "wall_sd": ok.std(ddof=1),
            "peak_mem_max": grouped["peak_mem_mb"].max(),
        }
    )
    summary["n_ok"] = summary["n_ok"].fillna(0).astype(int)
    summary["wall_cv"] = summary.pop("wall_sd") / summary["wall_mean"]
    columns = ["n", "n_ok", "wall_mean", "wall_cv", "peak_mem_max"]
    return summary[columns].reset_index()
