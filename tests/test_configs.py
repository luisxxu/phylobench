from pathlib import Path

import pytest
import yaml
from snakemake.utils import validate

from phylobench.emp import EMP_QUERY_N
from phylobench.regions import REGIONS

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
