"""The T28 overlays plan only column-window cells (pure planning, then dry runs)."""

from collections import Counter
from pathlib import Path

import pytest
import yaml
from snakemake.utils import validate
from test_workflow_dag import build_repo, planned_blocks

from phylobench.regions import COLUMN_REGIONS
from phylobench.tables import load_arms, plan_cells, run_specs

ROOT = Path(__file__).parents[1]
ARMS = load_arms(ROOT / "config" / "arms")
STAGES = {"S1": ("t28-b10k.yaml", "B10k"), "S2": ("t28-bfull.yaml", "Bfull")}
T28_ARMS = {"sepp-prod", "sepp-pinned", "bscampp"}


def merged(*names):
    config = {}
    for name in names:
        config.update(yaml.safe_load((ROOT / "config" / name).read_text()))
    return config


def cells(*names):
    return plan_cells(merged(*names), ARMS)


@pytest.mark.parametrize("name", ["t28.yaml", *(s[0] for s in STAGES.values())])
def test_t28_configs_validate(name):
    config = merged("t28.yaml", name)
    config["arms"] = ARMS
    validate(config, str(ROOT / "config" / "schema.yaml"))


def test_t28_plans_eight_column_regions_three_arms_two_backbones():
    planned = cells("t28.yaml")
    assert len(planned) == 48
    assert {c.region for c in planned} == set(COLUMN_REGIONS)
    assert {c.arm for c in planned} == T28_ARMS
    assert {c.aln for c in planned} == {"hmmalign", "raw"}
    assert {c.dataset for c in planned} == {"gg2"}


@pytest.mark.parametrize("stage", STAGES)
def test_each_stage_is_one_backbone_with_24_cells(stage):
    name, backbone = STAGES[stage]
    planned = cells("t28.yaml", name)
    assert Counter(c.backbone for c in planned) == {backbone: 24}


def test_stages_are_disjoint_and_cover_the_t28_matrix():
    s1, s2 = (set(cells("t28.yaml", s[0])) for s in STAGES.values())
    assert not s1 & s2
    assert s1 | s2 == set(cells("t28.yaml"))


def test_t28_never_plans_a_cell_t22_planned():
    assert not set(cells("t28.yaml")) & set(cells("full.yaml"))


@pytest.mark.parametrize("stage", STAGES)
def test_stage_run_specs_are_unique(stage):
    keys = [
        (s.run_id, s.phase)
        for s in run_specs(merged("t28.yaml", STAGES[stage][0]), ARMS)
    ]
    assert len(keys) == len(set(keys))


@pytest.fixture(scope="module")
def blocks_by_stage(tmp_path_factory):
    repo = tmp_path_factory.mktemp("t28")
    build_repo(repo)
    return {
        stage: planned_blocks(
            repo,
            "--configfile",
            "config/t28.yaml",
            f"config/{STAGES[stage][0]}",
        )
        for stage in STAGES
    }


@pytest.mark.parametrize("stage", STAGES)
def test_dry_run_plans_only_the_new_cells(blocks_by_stage, stage):
    blocks = blocks_by_stage[stage]
    backbone = STAGES[stage][1]
    other = next(b for _, b in STAGES.values() if b != backbone)
    assert len(blocks["extract_regions"]) == 8
    for rule in ("hmmalign", "sepp_prod", "sepp_pinned_place", "bscampp_place"):
        assert len(blocks[rule]) == 8, rule
    for rule in ("oracle_region", "pplacer_place", "epa_ng_place", "dipper_m_place"):
        assert rule not in blocks
    planned = "\n".join(b for rule in blocks.values() for b in rule)
    assert f"gg2/{other}/V" not in planned
    assert f"gg2/{backbone}/V" in planned


def test_dry_run_column_extraction_reads_the_shipped_alignment(blocks_by_stage):
    for block in blocks_by_stage["S1"]["extract_regions"]:
        assert "--aligned data/gg2/queries.aligned.fasta" in block
