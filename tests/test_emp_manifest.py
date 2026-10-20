import pandas as pd
import pytest
from emp_helpers import (
    BACKBONES,
    DATASETS,
    REGION,
    run_manifest,
    write_emp_stage,
)

from phylobench.similarity import BIN_LABELS, TOP_BIN


@pytest.fixture
def manifest(tmp_path):
    write_emp_stage(tmp_path)
    result = run_manifest(tmp_path)
    assert result.returncode == 0, result.stderr
    return pd.read_csv(
        tmp_path / "data" / "emp" / "MANIFEST.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


def test_columns_follow_t18_plus_dataset(manifest):
    assert list(manifest.columns) == (
        "dataset kind backbone region bin path n sha256 flagged".split()
    )


def test_file_rows_cover_biom_queries_and_backbone_links(manifest):
    files = manifest[manifest["kind"] == "file"]
    biom = files[files["path"].str.endswith(".biom")].iloc[0]
    assert biom["n"] == "3" and len(biom["sha256"]) == 64
    queries = files[files["path"].str.endswith("queries.fasta")]
    assert sorted(queries["dataset"]) == sorted(DATASETS)
    assert (queries["n"] == "2").all() and (queries["region"] == REGION).all()
    linked = files[files["backbone"] != ""]
    assert len(linked) == len(DATASETS) * len(BACKBONES) * 3
    assert (linked[~linked["path"].str.endswith("raxml-info.txt")]["n"] == "7").all()


def test_bin_rows_include_zero_bins_per_dataset_and_backbone(manifest):
    bins = manifest[manifest["kind"] == "bin"]
    assert len(bins) == len(DATASETS) * len(BACKBONES) * (len(BIN_LABELS) + 1)
    one = bins[(bins["dataset"] == "emp-10k") & (bins["backbone"] == "B10k")]
    assert one["bin"].tolist() == [*BIN_LABELS, TOP_BIN]
    assert one["n"].tolist() == list("000020000")
    assert (one["flagged"] == "1").all()
    assert (one["region"] == REGION).all()


def test_sha_mismatch_against_gg2_manifest_fails(tmp_path):
    write_emp_stage(tmp_path)
    (tmp_path / "data" / "gg2" / "B10k" / "tree.nwk").write_text("changed\n")
    result = run_manifest(tmp_path)
    assert result.returncode != 0
    assert "sha256 differs" in result.stderr


def test_backbone_file_missing_from_gg2_manifest_fails(tmp_path):
    write_emp_stage(tmp_path)
    manifest = tmp_path / "data" / "gg2" / "MANIFEST.tsv"
    table = pd.read_csv(manifest, sep="\t")
    table[~table["path"].str.endswith("raxml-info.txt")].to_csv(
        manifest, sep="\t", index=False
    )
    result = run_manifest(tmp_path)
    assert result.returncode != 0
    assert "no file row" in result.stderr
