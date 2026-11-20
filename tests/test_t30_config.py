"""config/t30.yaml plans the Bfull* stage on gg2-wol (planning, then a dry run)."""

from collections import Counter
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
CORE = {("raw", "sepp-prod"), ("raw", "sepp-pinned"), ("hmmalign", "bscampp")}
COST = {
    ("hmmalign", arm)
    for arm in (
        "pplacer-scampp",
        "epa-ng-scampp",
        "bscampp-b5000",
        "apples2",
        "dipper-m",
    )
}
HEADLINE_REGIONS = {"full-length", "V4-150"}
CELL_COUNTS = {
    "hmmalign": 9,
    "mask_hmmalign": 9,
    "similarity": 9,
    "sepp_prod": 9,
    "sepp_pinned_place": 9,
    "bscampp_place": 9,
    "pplacer_scampp_place": 2,
    "epa_ng_scampp_place": 2,
    "bscampp_b5000_place": 2,
    "apples2_place": 2,
    "dipper_m_place": 2,
    "dipper_m_input": 2,
    "dipper_to_jplace": 2,
    "guppy_tog": 26,
    "evaluate": 37,
}
ABSENT = (
    "upp_align",
    "upp_queries",
    "oracle_region",
    "pplacer_place",
    "epa_ng_place",
    "dipper_r_place",
    "dipper_m_gpu_place",
)


def load(name):
    return yaml.safe_load((ROOT / "config" / name).read_text())


def test_t30_config_validates():
    validate({**load("t30.yaml"), "arms": ARMS}, str(ROOT / "config" / "schema.yaml"))


def test_t30_plans_27_core_and_10_cost_cells_on_gg2_wol_bfull_only():
    cells = plan_cells(load("t30.yaml"), ARMS)
    assert len(cells) == 37
    assert {(c.dataset, c.backbone) for c in cells} == {("gg2-wol", "Bfull")}
    core = [c for c in cells if (c.aln, c.arm) in CORE]
    cost = [c for c in cells if (c.aln, c.arm) in COST]
    assert len(core) == 27 and len(cost) == 10
    assert Counter(c.arm for c in core) == {
        "sepp-prod": 9,
        "sepp-pinned": 9,
        "bscampp": 9,
    }
    assert {c.region for c in core} == REGIONS
    assert {c.region for c in cost} == HEADLINE_REGIONS
    assert {c.aln for c in cells} <= {"raw", "hmmalign"}


def test_t30_plans_no_excluded_arm():
    arms = {c.arm for c in plan_cells(load("t30.yaml"), ARMS)}
    assert not arms & {"pplacer", "epa-ng", "dipper-r", "dipper-m-gpu"}


def test_t30_run_ids_overlap_no_earlier_stage():
    t30 = {c.run_id for c in plan_cells(load("t30.yaml"), ARMS)}
    earlier = [load("full.yaml"), load("t28.yaml"), load("t29.yaml")]
    earlier.append(load("t29.yaml") | load("t29-combo.yaml"))
    for config in earlier:
        assert not t30 & {c.run_id for c in plan_cells(config, ARMS)}


def test_run_specs_are_unique_with_one_hmmbuild_and_nine_hmmaligns():
    specs = run_specs(load("t30.yaml"), ARMS)
    keys = [(s.run_id, s.phase) for s in specs]
    assert len(keys) == len(set(keys))
    aligns = [s for s in specs if s.phase == "align"]
    assert len([s for s in aligns if s.aln == "hmmalign"]) == 9
    assert len([s for s in aligns if s.region == ""]) == 1
    assert len([s for s in specs if s.phase == "ref-prep"]) == 1


@pytest.fixture(scope="module")
def blocks(tmp_path_factory):
    repo = tmp_path_factory.mktemp("t30")
    build_repo(repo)
    return planned_blocks(repo, "--configfile", "config/t30.yaml")


def test_dry_run_cell_counts(blocks):
    assert {rule: len(blocks[rule]) for rule in CELL_COUNTS} == CELL_COUNTS


def test_dry_run_plans_no_upp_oracle_or_excluded_arm(blocks):
    assert not [rule for rule in ABSENT if rule in blocks]


def test_dry_run_place_blocks_read_the_gg2_wol_bfull_backbone(blocks):
    for rule in (
        "sepp_prod",
        "sepp_pinned_place",
        "bscampp_place",
        "pplacer_scampp_place",
        "epa_ng_scampp_place",
        "bscampp_b5000_place",
        "apples2_place",
        "dipper_m_place",
    ):
        assert all("data/gg2-wol/Bfull/" in block for block in blocks[rule]), rule
