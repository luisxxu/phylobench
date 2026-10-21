from pathlib import Path

import pytest
import yaml
from snakemake.utils import validate

from phylobench.emp import EMP_QUERY_N
from phylobench.regions import REGIONS
from phylobench.tables import plan_cells

ROOT = Path(__file__).parents[1]
SCHEMA = ROOT / "config" / "schema.yaml"
ARM_FILES = sorted((ROOT / "config" / "arms").glob("*.yaml"))


def load(path):
    return yaml.safe_load(path.read_text())


def test_schema_regions_match_regions_module():
    enum = load(SCHEMA)["properties"]["regions"]["items"]["enum"]
    assert tuple(enum) == REGIONS


def test_smoke_config_validates_with_every_arm():
    config = load(ROOT / "config" / "smoke.yaml")
    config["arms"] = {path.stem: load(path) for path in ARM_FILES}
    validate(config, str(SCHEMA))


@pytest.mark.parametrize("path", ARM_FILES, ids=lambda p: p.stem)
def test_arm_files_validate(path):
    validate({"arms": {path.stem: load(path)}}, str(SCHEMA))


@pytest.mark.parametrize("missing", ["tool", "aligned", "phases"])
def test_arm_requires_tool_aligned_phases(missing):
    arm = {"tool": "x", "aligned": True, "phases": {"place": "place.run.json"}}
    del arm[missing]
    with pytest.raises(Exception, match=missing):
        validate({"arms": {"x": arm}}, str(SCHEMA))


def test_arm_rejects_unknown_keys():
    arm = {
        "tool": "x",
        "aligned": True,
        "phases": {"p": "p.json"},
        "skip_backbones": [],
    }
    with pytest.raises(Exception, match="skip_backbones"):
        validate({"arms": {"x": arm}}, str(SCHEMA))


def test_emp_config_validates_and_is_emp_only():
    config = load(ROOT / "config" / "emp.yaml")
    config["arms"] = {path.stem: load(path) for path in ARM_FILES}
    validate(config, str(SCHEMA))
    assert set(config["datasets"]) == set(EMP_QUERY_N)
    assert config["regions"] == ["V4-150"] and config["alns"] == ["hmmalign"]


def pilot_config(name):
    config = load(ROOT / "config" / name)
    config["arms"] = {path.stem: load(path) for path in ARM_FILES}
    return config


@pytest.mark.parametrize(
    "name", ["pilot.yaml", "pilot-cv.yaml", "pilot-bfull-probe.yaml"]
)
def test_pilot_configs_validate(name):
    validate(pilot_config(name), str(SCHEMA))


def test_pilot_plans_every_arm_in_33_cells():
    config = pilot_config("pilot.yaml")
    cells = plan_cells(config, config["arms"])
    assert len(cells) == 33
    assert {c.arm for c in cells} == {path.stem for path in ARM_FILES}


def test_pilot_cv_plans_one_cell_per_arm():
    config = pilot_config("pilot-cv.yaml")
    assert len(plan_cells(config, config["arms"])) == len(ARM_FILES) == 11


def test_timeout_defaults_to_71_hours():
    config = {"arms": {}}
    validate(config, str(SCHEMA))
    assert config["timeout_s"] == 255600
