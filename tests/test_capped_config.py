from pathlib import Path

import yaml
from snakemake.utils import validate

from phylobench.tables import load_arms, plan_cells

ROOT = Path(__file__).parents[1]
SCHEMA = ROOT / "config" / "schema.yaml"


def load(path):
    return yaml.safe_load(path.read_text())


def capped_config():
    config = load(ROOT / "config" / "smoke.yaml") | load(
        ROOT / "config" / "pplacer-bfull-cap.yaml"
    )
    config["arms"] = load_arms(ROOT / "config" / "arms", only=config.pop("arm"))
    return config


def test_overlay_validates_against_schema():
    validate(capped_config(), str(SCHEMA))


def test_overlay_plans_exactly_one_pplacer_cell():
    config = capped_config()
    cells = plan_cells(config, config["arms"])
    assert [(c.backbone, c.region, c.aln, c.arm) for c in cells] == [
        ("Bfull", "V4-150", "hmmalign", "pplacer")
    ]
