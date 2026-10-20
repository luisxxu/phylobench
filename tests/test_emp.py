import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from phylobench.emp import EMP_QUERY_N, asv_failures, asv_names, query_set

ROOT = Path(__file__).parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
GOOD = "ACGT" * 37 + "AC"


def make_ids(count):
    rng = np.random.default_rng(1)
    bases = np.array(list("ACGT"))
    return ["".join(bases[rng.integers(4, size=150)]) for _ in range(count)]


def test_good_ids_have_no_failures():
    assert asv_failures([GOOD, GOOD.replace("A", "C")]) == []


@pytest.mark.parametrize(
    "bad",
    [
        GOOD[:-1] + "N",
        GOOD[:-1] + "U",
        GOOD[:-1] + "a",
        GOOD[:-1],
        GOOD + "A",
        GOOD + "\n",
        "",
    ],
)
def test_bad_ids_are_reported(bad):
    assert asv_failures([GOOD, bad]) == [bad]


def test_duplicates_are_reported_once():
    assert asv_failures([GOOD, GOOD, GOOD]) == [GOOD]


def test_query_set_is_deterministic_and_order_independent():
    ids = make_ids(200)
    assert query_set(ids, 50) == query_set(ids, 50)
    assert query_set(ids, 50) == query_set(ids[::-1], 50)
    assert query_set(ids, 50, seed=1) != query_set(ids, 50, seed=0)


def test_query_sets_are_nested_unique_subsets():
    ids = make_ids(200)
    large, small = query_set(ids, 100), query_set(ids, 25)
    assert small == large[:25]
    assert len(set(large)) == 100 and set(large) <= set(ids)


def test_query_set_rejects_n_above_available():
    with pytest.raises(ValueError, match="only 3"):
        query_set(["a", "b", "c"], 4)


def test_asv_names_are_short_and_ordered():
    assert asv_names(3) == ["asv000001", "asv000002", "asv000003"]
    assert asv_names(EMP_QUERY_N["emp-100k"])[-1] == "asv100000"


def write_biom(path, ids, version=(2, 1)):
    with h5py.File(path, "w") as f:
        f.attrs["format-version"] = np.array(version)
        f.create_dataset("observation/ids", data=np.array(ids, dtype="S"))


def run_query_sets(tmp_path, ids, leaves, n=5, version=(2, 1)):
    biom, tree = tmp_path / "x.biom", tmp_path / "tree.nwk"
    write_biom(biom, ids, version)
    tree.write_text("(" + ",".join(f"{leaf}:1" for leaf in leaves) + ");")
    return subprocess.run(
        [
            sys.executable, str(SCRIPTS / "emp_query_sets.py"),
            "--biom", str(biom), "--backbone-tree", str(tree),
            "--n", str(n), "--out", str(tmp_path / "queries.fasta"),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )  # fmt: skip


def test_script_writes_named_fasta_of_sampled_asvs(tmp_path):
    ids = make_ids(20)
    result = run_query_sets(tmp_path, ids, ["a", "b", "c"])
    assert result.returncode == 0, result.stderr
    lines = (tmp_path / "queries.fasta").read_text().splitlines()
    assert lines[0::2] == [f">{name}" for name in asv_names(5)]
    assert lines[1::2] == query_set(ids, 5)


def test_script_rejects_unexpected_format_version(tmp_path):
    result = run_query_sets(tmp_path, make_ids(20), ["a"], version=(1, 0))
    assert result.returncode != 0
    assert "format-version" in result.stderr
    assert not (tmp_path / "queries.fasta").exists()


def test_script_rejects_non_sequence_ids(tmp_path):
    result = run_query_sets(tmp_path, [*make_ids(20), "nonsense"], ["a"])
    assert result.returncode != 0
    assert "nonsense" in result.stderr


def test_script_rejects_name_colliding_with_backbone_leaf(tmp_path):
    result = run_query_sets(tmp_path, make_ids(20), ["asv000003", "b"])
    assert result.returncode != 0
    assert "asv000003" in result.stderr
