"""config/t28.yaml plans only gg2-wol B10k* cells (pure planning, then a dry run)."""

from pathlib import Path

import pytest
import yaml
from snakemake.utils import validate
from test_workflow_dag import build_repo, planned_blocks

from phylobench.tables import load_arms, plan_cells, run_specs

ROOT = Path(__file__).parents[1]
ARMS = load_arms(ROOT / "config" / "arms")
REGIONS = {
    "full-length",
    "V4-150",
    "V3-V4",
    "V1-V2",
    "V1-V3",
    "V4-V5",
    "V5-V6",
    "V6-V8",
    "V7-V9c",
}
T28_ARMS = {"sepp-prod", "sepp-pinned", "bscampp", "pplacer-scampp", "epa-ng-scampp"}


def load(name):
    return yaml.safe_load((ROOT / "config" / name).read_text())


def test_t28_config_validates():
    config = {**load("t28.yaml"), "arms": ARMS}
    validate(config, str(ROOT / "config" / "schema.yaml"))


def test_t28_plans_nine_regions_five_arms_at_one_backbone():
    planned = plan_cells(load("t28.yaml"), ARMS)
    assert len(planned) == 45
    assert {c.region for c in planned} == REGIONS
    assert {c.arm for c in planned} == T28_ARMS
    assert {(c.dataset, c.backbone) for c in planned} == {("gg2-wol", "B10k")}
    assert {c.aln for c in planned} == {"hmmalign", "raw"}


def test_t28_never_plans_a_cell_t22_planned():
    t22 = set(plan_cells(load("full.yaml"), ARMS))
    assert not set(plan_cells(load("t28.yaml"), ARMS)) & t22


def test_run_specs_are_unique():
    keys = [(s.run_id, s.phase) for s in run_specs(load("t28.yaml"), ARMS)]
    assert len(keys) == len(set(keys))


@pytest.fixture(scope="module")
def blocks(tmp_path_factory):
    repo = tmp_path_factory.mktemp("t28")
    build_repo(repo)
    return planned_blocks(repo, "--configfile", "config/t28.yaml")


def test_dry_run_plans_each_region_once_and_each_cell_once(blocks):
    assert len(blocks["extract_regions"]) == 9
    assert len(blocks["hmmalign"]) == 9
    for rule in (
        "sepp_prod",
        "sepp_pinned_place",
        "bscampp_place",
        "pplacer_scampp_place",
        "epa_ng_scampp_place",
    ):
        assert len(blocks[rule]) == 9, rule
    assert len(blocks["evaluate"]) == 45
    for rule in ("oracle_region", "pplacer_place", "epa_ng_place", "dipper_m_place"):
        assert rule not in blocks


def test_dry_run_never_touches_gg2_or_bfull_cells(blocks):
    cells = "\n".join(
        b for rule, bs in blocks.items() if rule != "gg2_wol_build" for b in bs
    )
    assert "results/gg2-wol/B10k/" in cells
    assert "results/gg2/" not in cells
    assert "Bfull" not in cells.replace("data/gg2-wol/Bfull", "")


def test_dry_run_extracts_the_column_region_from_the_wol_alignment(blocks):
    (block,) = [b for b in blocks["extract_regions"] if "region V7-V9c" in b]
    assert "--aligned data/gg2-wol/queries.aligned.fasta" in block
    others = [b for b in blocks["extract_regions"] if "region V7-V9c" not in b]
    assert not any("--aligned" in b for b in others)
