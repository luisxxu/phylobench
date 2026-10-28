from pathlib import Path

import pandas as pd
import pytest

from phylobench.tables import ACCURACY_COLUMNS, RUNS_COLUMNS

SAMPLE = Path(__file__).parent / "fixtures" / "sample_results"


@pytest.fixture(scope="module")
def runs():
    return pd.read_csv(SAMPLE / "runs.tsv", sep="\t")


@pytest.fixture(scope="module")
def accuracy():
    return pd.read_csv(SAMPLE / "accuracy.tsv", sep="\t")


def test_columns_match_the_table_constants(runs, accuracy):
    assert tuple(runs.columns) == RUNS_COLUMNS
    assert tuple(accuracy.columns) == ACCURACY_COLUMNS


def test_files_stay_small():
    assert all(
        (SAMPLE / f).stat().st_size < 1_000_000 for f in ("runs.tsv", "accuracy.tsv")
    )


def test_every_accuracy_run_has_an_evaluate_row(runs, accuracy):
    evaluated = set(runs.loc[runs["phase"] == "evaluate", "run_id"])
    assert set(accuracy["run_id"]) == evaluated


def test_covers_datasets_backbones_regions_and_all_arms(runs):
    assert {"gg2", "sim10k", "emp-10k", "emp-100k"} <= set(runs["dataset"])
    assert {"B10k", "Bfull"} == set(runs["backbone"])
    assert {"full-length", "V4-150"} <= set(runs["region"])
    assert runs["arm"].nunique() == 11
    assert runs.loc[runs["aln"] == "oracle", "backbone"].eq("B10k").all()


def test_edge_cases_are_present(runs, accuracy):
    dnf = runs[runs["status"] == "DNF"]
    apples = dnf[(dnf["arm"] == "apples2") & (dnf["phase"] == "convert")]
    assert len(apples) >= 2
    assert (dnf["dnf_reason"] == "timeout").any()
    assert (runs["mem_source"] == "psutil").any()
    assert runs["cpu_model"].str.contains("Xeon").any()
    assert runs.loc[runs["arm"] == "dipper-m-gpu", "gpu_peak_mb"].notna().any()
    assert (accuracy["bin"] == "<75").any()
    pplacer = accuracy[accuracy["run_id"] == "gg2/Bfull/V4-150/hmmalign/pplacer"]
    assert len(pplacer) and pplacer["edge_dist"].isna().all()
    bfull = accuracy[accuracy["run_id"].str.contains("/Bfull/")]
    assert bfull["delta_err"].isna().any()
