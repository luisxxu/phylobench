"""Run planning and the tidy runs/accuracy tables (docs/harness.md section 7)."""

import math
from collections import Counter
from itertools import product
from pathlib import Path
from typing import NamedTuple

import yaml
from snakemake.exceptions import WorkflowError

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
FORCED_DNF = "forced-dnf"
RAW = "raw"
FULL_LENGTH = "full-length"


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


def plan_cells(config: dict, arms: dict[str, dict]) -> list[Cell]:
    cells = []
    for dataset, backbone, region, (name, arm) in product(
        config["datasets"], config["backbones"], config["regions"], arms.items()
    ):
        for aln in config["alns"] if arm["aligned"] else [RAW]:
            cells.append(Cell(dataset, backbone, region, aln, name))
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
    specs = []
    for dataset, backbone in product(config["datasets"], config["backbones"]):
        base = f"results/{dataset}/{backbone}"
        shared = f"{dataset}/{backbone}"
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
        for region in config["regions"]:
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
    for cell in plan_cells(config, arms):
        arm = arms[cell.arm]
        for phase, filename in arm["phases"].items():
            tool = "guppy" if phase == "convert" else arm["tool"]
            specs.append(cell_spec(cell, f"{cell.dir}/{filename}", phase, tool))
        specs.append(
            cell_spec(cell, f"{cell.dir}/evaluate.run.json", "evaluate", "phylobench")
        )
    for cell in forced_dnf_cell(config):
        specs.append(cell_spec(cell, f"{cell.dir}/place.run.json", "place", "sleep"))
    return specs


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


def smoke_failures(runs, accuracy, config: dict, arms: dict, queries: list[str]):
    failures = []
    specs = run_specs(config, arms)
    planned = {(s.run_id, s.phase) for s in specs}
    seen = Counter(zip(runs["run_id"], runs["phase"], strict=True))
    for key in sorted(planned):
        if seen[key] != 1:
            failures.append(f"{key}: expected exactly one row, found {seen[key]}")
    for key in sorted(set(seen) - planned):
        failures.append(f"{key}: unplanned row")

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
        if sorted(scored["query"]) != wanted:
            failures.append(
                f"{cell.run_id}: {len(scored)} accuracy rows, queries differ from "
                f"the {len(wanted)} fixture queries"
            )
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
