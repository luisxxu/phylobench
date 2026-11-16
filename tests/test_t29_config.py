"""config/t29.yaml plans the B10k* lever stage on gg2-wol (planning, then a dry run)."""

from pathlib import Path

import pytest
import yaml
from snakemake.utils import validate
from test_workflow_dag import build_repo, planned_blocks

from phylobench.tables import load_arms, plan_cells, run_specs

ROOT = Path(__file__).parents[1]
ARMS = load_arms(ROOT / "config" / "arms")
REUSED = {("raw", "sepp-prod"), ("hmmalign", "bscampp")}
NEW = {("upp", "bscampp"), ("hmmalign", "bscampp-b5000")}


def load(*names):
    config = {}
    for name in names:
        config.update(yaml.safe_load((ROOT / "config" / name).read_text()))
    return config


def kinds(cells):
    return {(c.aln, c.arm) for c in cells}


def test_t29_configs_validate():
    for names in (("t29.yaml",), ("t29.yaml", "t29-combo.yaml")):
        validate({**load(*names), "arms": ARMS}, str(ROOT / "config" / "schema.yaml"))


def test_t29_plans_36_cells_and_the_combo_adds_nine():
    cells = plan_cells(load("t29.yaml"), ARMS)
    assert len(cells) == 36
    assert kinds(cells) == REUSED | NEW
    assert {(c.dataset, c.backbone) for c in cells} == {("gg2-wol", "B10k")}
    assert len({c.region for c in cells}) == 9
    combo = plan_cells(load("t29.yaml", "t29-combo.yaml"), ARMS)
    assert len(combo) == 45
    assert kinds(combo) == REUSED | NEW | {("upp", "bscampp-b5000")}


def test_reused_cells_are_t28_cells_with_identical_run_ids():
    t28 = {c.run_id for c in plan_cells(load("t28.yaml"), ARMS)}
    t29 = plan_cells(load("t29.yaml"), ARMS)
    reused = {c.run_id for c in t29 if (c.aln, c.arm) in REUSED}
    assert len(reused) == 18 and reused <= t28
    assert not {c.run_id for c in t29 if (c.aln, c.arm) in NEW} & t28


def test_t29_never_plans_a_t22_cell_or_a_bfull_cell():
    t22 = set(plan_cells(load("full.yaml"), ARMS))
    cells = plan_cells(load("t29.yaml", "t29-combo.yaml"), ARMS)
    assert not set(cells) & t22
    assert not [c for c in cells if c.dataset == "gg2" or c.backbone == "Bfull"]


def test_run_specs_are_unique_and_include_each_upp_region_once():
    specs = run_specs(load("t29.yaml"), ARMS)
    keys = [(s.run_id, s.phase) for s in specs]
    assert len(keys) == len(set(keys))
    assert len([s for s in specs if s.aln == "upp" and s.phase == "align"]) == 9


@pytest.fixture(scope="module")
def blocks(tmp_path_factory):
    repo = tmp_path_factory.mktemp("t29")
    build_repo(repo)
    return planned_blocks(repo, "--configfile", "config/t29.yaml")


def test_dry_run_plans_the_upp_chain_once_per_region(blocks):
    assert len(blocks["upp_align"]) == 9
    assert len(blocks["upp_queries"]) == 9
    for block in blocks["upp_align"]:
        assert "-a data/gg2-wol/B10k/aligned-dna-sequences.fasta" in block
        assert "-t data/gg2-wol/B10k/tree.nwk" in block
        assert "-A 100" in block
        assert block.rindex("--timeout 84000") > block.index("--timeout 255600")


def test_dry_run_places_each_new_cell_once_from_its_own_alignment(blocks):
    placed = blocks["bscampp_place"]
    assert len(placed) == 18
    assert sum("upp.aligned.fasta" in b for b in placed) == 9
    b5000 = blocks["bscampp_b5000_place"]
    assert len(b5000) == 9
    assert all("-b 5000 -V 5" in b for b in b5000)
    assert all("hmmalign.aligned.fasta" in b for b in b5000)
    # nothing exists in the test repo, so the 18 reused cells are planned too
    assert len(blocks["guppy_tog"]) == 27 and len(blocks["evaluate"]) == 36


def test_dry_run_plans_no_other_arm(blocks):
    for rule in (
        "pplacer_scampp_place",
        "epa_ng_scampp_place",
        "sepp_pinned_place",
        "oracle_region",
    ):
        assert rule not in blocks


def test_combo_adds_b5000_on_upp(tmp_path):
    build_repo(tmp_path)
    combo = planned_blocks(
        tmp_path, "--configfile", "config/t29.yaml", "config/t29-combo.yaml"
    )
    b5000 = combo["bscampp_b5000_place"]
    assert len(b5000) == 18
    assert sum("upp.aligned.fasta" in b for b in b5000) == 9
