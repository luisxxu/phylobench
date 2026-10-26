"""The T22 matrix, planned in pure Python (no dry runs)."""

from collections import Counter
from pathlib import Path

import pytest
import yaml
from snakemake.utils import validate

from phylobench.emp import EMP_QUERY_N, REFERENCE_ARM
from phylobench.tables import load_arms, plan_cells, run_specs

ROOT = Path(__file__).parents[1]
SCHEMA = ROOT / "config" / "schema.yaml"
ARMS = load_arms(ROOT / "config" / "arms")
STAGES = {"S1": "full-b10k.yaml", "S2": "full-bfull.yaml"}


def merged(*names):
    config = {}
    for name in names:
        config.update(yaml.safe_load((ROOT / "config" / name).read_text()))
    return config


def cells(*names, arms=ARMS):
    return plan_cells(merged(*names), arms)


@pytest.mark.parametrize("name", ["full.yaml", *STAGES.values()])
def test_full_configs_validate(name):
    config = merged("full.yaml", name)
    config["arms"] = ARMS
    validate(config, str(SCHEMA))


def test_counts_per_dataset_and_backbone():
    counts = Counter((c.dataset, c.backbone) for c in cells("full.yaml"))
    assert counts == {
        ("gg2", "B10k"): 171,
        ("sim10k", "B10k"): 22,
        ("emp-10k", "B10k"): 11,
        ("emp-100k", "B10k"): 11,
        ("gg2", "Bfull"): 73,
        ("emp-10k", "Bfull"): 8,
        ("emp-100k", "Bfull"): 8,
    }
    assert len(cells("full.yaml")) == 304


def test_stages_are_disjoint_and_cover_the_full_matrix():
    s1 = set(cells("full.yaml", STAGES["S1"]))
    s2 = set(cells("full.yaml", STAGES["S2"]))
    assert (len(s1), len(s2)) == (215, 89)
    assert not s1 & s2
    assert s1 | s2 == set(cells("full.yaml"))


def test_bfull_runs_only_what_is_tractable():
    bfull = [c for c in cells("full.yaml") if c.backbone == "Bfull"]
    assert not [c for c in bfull if c.arm in ("epa-ng", "dipper-r")]
    assert not [c for c in bfull if c.aln == "oracle"]
    assert not [c for c in bfull if c.dataset == "sim10k"]
    pplacer = [c for c in bfull if c.arm == "pplacer"]
    assert [c.run_id for c in pplacer] == ["gg2/Bfull/V4-150/hmmalign/pplacer"]


def test_capped_pplacer_overlay_is_the_one_bfull_pplacer_cell():
    capped = cells(
        "full.yaml", "pplacer-bfull-cap.yaml", arms={"pplacer": ARMS["pplacer"]}
    )
    assert [c.run_id for c in capped] == ["gg2/Bfull/V4-150/hmmalign/pplacer"]
    assert capped[0] in cells("full.yaml")


def test_sim10k_runs_two_regions_hmmalign_all_arms_at_b10k():
    sim = [c for c in cells("full.yaml") if c.dataset == "sim10k"]
    assert {(c.backbone, c.region) for c in sim} == {
        ("B10k", "full-length"),
        ("B10k", "V4-150"),
    }
    assert {c.aln for c in sim} == {"hmmalign", "raw"}
    assert {c.arm for c in sim} == set(ARMS)


@pytest.mark.parametrize("stage", [None, *STAGES])
def test_dipper_m_and_gpu_plan_identical_cells(stage):
    names = ["full.yaml", *([STAGES[stage]] if stage else [])]

    def cell_set(arm):
        return {c[:4] for c in cells(*names) if c.arm == arm}

    assert cell_set("dipper-m") == cell_set("dipper-m-gpu")


def test_sepp_prod_covers_every_emp_dataset_and_backbone():
    config = merged("full.yaml")
    planned = cells("full.yaml")
    for dataset in EMP_QUERY_N:
        for backbone in config["backbones"]:
            assert any(
                (c.dataset, c.backbone, c.arm) == (dataset, backbone, REFERENCE_ARM)
                for c in planned
            )


def test_every_stage_run_spec_is_unique():
    for name in STAGES.values():
        config = merged("full.yaml", name)
        keys = [(s.run_id, s.phase) for s in run_specs(config, ARMS)]
        assert len(keys) == len(set(keys))
