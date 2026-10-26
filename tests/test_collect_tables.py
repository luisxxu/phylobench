import json

import pandas as pd
from emp_helpers import run_script, write

from phylobench.tables import (
    ACCURACY_COLUMNS,
    CONCORDANCE_COLUMNS,
    RunSpec,
    concordance_rows,
)

RUN = {
    "status": "ok",
    "wall_s": 1.0,
    "cpu_s": 1.0,
    "peak_mem_mb": 1.0,
    "mem_source": "psutil",
    "threads": 1,
    "node": "n",
    "cpu_model": "c",
    "version": "v",
    "env_lock_sha256": "e",
}
GG2_ID = "gg2/B10k/V4-150/hmmalign/pplacer"
EMP_ID = "emp-10k/B10k/V4-150/hmmalign/pplacer"


def spec(run_id, json_name):
    dataset, backbone, region, aln, arm = run_id.split("/")
    return RunSpec(
        f"results/{run_id}/{json_name}",
        run_id,
        "evaluate",
        "phylobench",
        dataset,
        backbone,
        region,
        aln,
        arm,
    )


def stage(tmp_path):
    """Both datasets in one stage: a gg2 accuracy cell and an EMP concordance cell."""
    specs = [spec(GG2_ID, "evaluate.run.json"), spec(EMP_ID, "concordance.run.json")]
    for s in specs:
        write(tmp_path / s.run_json, json.dumps(RUN))
        write(tmp_path / f"results/{s.dataset}/regions/V4-150/queries.fasta", ">q\nA\n")
    pd.DataFrame(
        [
            {
                "run_id": GG2_ID,
                "query": "q",
                "sim": 0.9,
                "bin": "90-95",
                "edge_dist": 0,
                "exact": 1,
                "delta_err": 0,
                "label": "accuracy",
            }
        ],
        columns=list(ACCURACY_COLUMNS),
    ).to_csv(tmp_path / f"results/{GG2_ID}/accuracy.tsv", sep="\t", index=False)
    sims = pd.DataFrame({"query": ["q"], "sim": [0.9], "bin": ["90-95"]})
    rows = concordance_rows(EMP_ID, "ref", [{"query": "q", "edge_dist": 2.0}], sims)
    pd.DataFrame(rows, columns=list(CONCORDANCE_COLUMNS)).to_csv(
        tmp_path / f"results/{EMP_ID}/concordance.tsv", sep="\t", index=False
    )
    return specs


def collect(tmp_path, specs):
    specs_file = write(
        tmp_path / "results/run_specs.json", json.dumps([s._asdict() for s in specs])
    )
    return run_script(
        "collect_tables.py",
        "--specs", specs_file,
        "--out-runs", tmp_path / "runs.tsv",
        "--out-accuracy", tmp_path / "accuracy.tsv",
        "--out-concordance", tmp_path / "concordance.tsv",
        cwd=tmp_path,
    )  # fmt: skip


def test_scores_are_split_by_table_kind(tmp_path):
    result = collect(tmp_path, stage(tmp_path))
    assert result.returncode == 0, result.stderr
    accuracy = pd.read_csv(tmp_path / "accuracy.tsv", sep="\t")
    concordance = pd.read_csv(tmp_path / "concordance.tsv", sep="\t")
    runs = pd.read_csv(tmp_path / "runs.tsv", sep="\t")
    assert accuracy["run_id"].tolist() == [GG2_ID]
    assert concordance["run_id"].tolist() == [EMP_ID]
    assert tuple(accuracy.columns) == ACCURACY_COLUMNS
    assert tuple(concordance.columns) == CONCORDANCE_COLUMNS
    assert sorted(runs["run_id"]) == sorted([GG2_ID, EMP_ID])


def test_no_emp_cells_gives_header_only_concordance(tmp_path):
    result = collect(tmp_path, stage(tmp_path)[:1])
    assert result.returncode == 0, result.stderr
    table = pd.read_csv(tmp_path / "concordance.tsv", sep="\t")
    assert table.empty and tuple(table.columns) == CONCORDANCE_COLUMNS
