"""Acceptance (T19): no output column labels concordance as accuracy."""

import re

import pandas as pd
import pytest
from emp_helpers import run_manifest, write_emp_stage
from test_collect_tables import collect, stage
from test_concordance_run import concordance

FORBIDDEN = re.compile("accuracy|truth|error", re.IGNORECASE)


def assert_not_accuracy(path):
    table = pd.read_csv(path, sep="\t")
    assert not [c for c in table.columns if FORBIDDEN.search(c)], path
    if "label" in table:
        assert set(table["label"]) <= {"concordance"}


@pytest.mark.parametrize("gate", ["ok", "timeout"])
def test_concordance_run_outputs(tmp_path, gate):
    assert concordance(tmp_path, gate=gate).returncode == 0
    assert_not_accuracy(tmp_path / "concordance.tsv")


def test_emp_manifest_columns(tmp_path):
    write_emp_stage(tmp_path)
    assert run_manifest(tmp_path).returncode == 0
    assert_not_accuracy(tmp_path / "data" / "emp" / "MANIFEST.tsv")


def test_collected_concordance_table(tmp_path):
    assert collect(tmp_path, stage(tmp_path)).returncode == 0
    table = pd.read_csv(tmp_path / "concordance.tsv", sep="\t")
    assert len(table) and set(table["label"]) == {"concordance"}
    assert_not_accuracy(tmp_path / "concordance.tsv")
    assert "concordance" not in set(
        pd.read_csv(tmp_path / "accuracy.tsv", sep="\t")["label"]
    )
