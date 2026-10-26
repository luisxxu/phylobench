"""Dry-run DAG checks. Only raw downloads exist; the gg2 rules are planned for real."""

import os
import re
import subprocess
import sys
from collections import defaultdict

import pytest
from emp_helpers import ROOT

LINKED = ("workflow", "envs", "config", "src", "tests", "profiles")
MIXED = """\
datasets: [gg2, emp-10k]
backbones: [B10k, Bfull]
regions: [V4-150]
alns: [hmmalign]
"""


def build_repo(tmp_path):
    for name in LINKED:
        (tmp_path / name).symlink_to(ROOT / name)
    for name in (
        "2022.10.backbone.sepp-reference.qza",
        "2022.10.backbone.tax.qza",
    ):
        (tmp_path / "data" / "gg2" / "download" / name).parent.mkdir(
            parents=True, exist_ok=True
        )
        (tmp_path / "data" / "gg2" / "download" / name).touch()
    biom = tmp_path / "data" / "emp" / "emp_deblur_150bp.release1.biom"
    biom.parent.mkdir(parents=True)
    biom.touch()
    (tmp_path / "mixed.yaml").write_text(MIXED)


@pytest.fixture
def repo(tmp_path):
    build_repo(tmp_path)
    return tmp_path


def planned_blocks(repo, *args) -> dict[str, list[str]]:
    result = subprocess.run(
        [
            sys.executable, "-m", "snakemake", "all", "-n", "-p",
            "--snakefile", "workflow/Snakefile", *args,
        ],
        capture_output=True,
        text=True,
        cwd=repo,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )  # fmt: skip
    assert result.returncode == 0, result.stderr[-3000:]
    blocks = defaultdict(list)
    for block in re.split(r"\n(?=(?:local)?rule \w+:)", result.stdout)[1:]:
        blocks[re.match(r"(?:local)?rule (\w+):", block)[1]].append(block)
    return blocks


def dry_run(repo, configfile) -> defaultdict[str, list[str]]:
    result = subprocess.run(
        [
            sys.executable, "-m", "snakemake", "-n", "-p",
            "--snakefile", "workflow/Snakefile", "--configfile", configfile,
        ],
        capture_output=True,
        text=True,
        cwd=repo,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )  # fmt: skip
    assert result.returncode == 0, result.stderr[-3000:]
    outputs = defaultdict(list)
    rule = None
    for line in result.stdout.splitlines():
        if match := re.fullmatch(r"(?:localrule|rule) (\w+):", line):
            rule = match[1]
        elif rule and (match := re.match(r"\s+output: (.*)", line)):
            outputs[rule].extend(match[1].split(", "))
    return defaultdict(list, {r: sorted(p) for r, p in outputs.items()})


@pytest.mark.parametrize(
    "configfile", ["config/smoke.yaml", "config/emp.yaml", "mixed.yaml"]
)
def test_dry_run_is_stable(repo, configfile):
    first = dry_run(repo, configfile)
    assert first == dry_run(repo, configfile)


def test_smoke_never_touches_emp(repo):
    planned = dry_run(repo, "config/smoke.yaml")
    assert not [r for r in planned if r.startswith("emp_") or r == "concordance"]
    assert (
        planned["evaluate"] and "results/concordance.tsv" in planned["collect_tables"]
    )


def test_emp_config_plans_concordance_not_evaluate(repo):
    planned = dry_run(repo, "config/emp.yaml")
    assert "evaluate" not in planned
    arms = len(list((ROOT / "config" / "arms").glob("*.yaml")))
    assert len(planned["concordance"]) == 2 * 2 * 2 * arms  # run.json + tsv per cell
    assert all("/emp-" in p for p in planned["concordance"])
    assert "data/emp/MANIFEST.tsv" in planned["emp_manifest"]
    assert "emp_download" not in planned and "gg2_download" not in planned
    assert planned["gg2_build"] and "data/gg2/MANIFEST.tsv" in planned["gg2_build"]
    assert len(planned["emp_query_sets"]) == 2
    assert len(planned["emp_hmm"]) == 4 and len(planned["emp_refprep"]) == 8
    assert not [p for p in planned["refprep_me_tree"] if "/emp-" in p]
    assert not [p for p in planned["hmmbuild"] if "/emp-" in p]
    assert not [p for p in planned["extract_regions"] if "/emp-" in p]


def test_mixed_config_keeps_evaluate_for_gg2_only(repo):
    planned = dry_run(repo, "mixed.yaml")
    assert planned["evaluate"]
    assert not [p for p in planned["evaluate"] if "/emp-" in p]
    assert not [p for p in planned["concordance"] if "/gg2/" in p]


def test_measure_commands_carry_global_and_per_rule_timeouts(repo):
    blocks = planned_blocks(
        repo, "-c", "16", "--configfile", "config/smoke.yaml", "--config", "threads=16"
    )
    assert "-m phylobench.measure --timeout 255600 --out" in blocks["sepp_prod"][0]
    pplacer = blocks["pplacer_place"][0]
    assert pplacer.rindex("--timeout 86400") > pplacer.index("--timeout 255600")


def test_gpu_arm_forces_binary_and_cpu_dipper_does_not(repo):
    blocks = planned_blocks(
        repo, "-c", "16", "--configfile", "config/smoke.yaml", "--config", "threads=16"
    )
    gpu, cpu = blocks["dipper_m_gpu_place"][0], blocks["dipper_m_place"][0]
    assert "--force-binary" in gpu and " -T " not in gpu
    assert "--force-binary" not in cpu and "-T 16" in cpu


@pytest.fixture(scope="module")
def full_blocks(tmp_path_factory):
    repo_dir = tmp_path_factory.mktemp("full")
    build_repo(repo_dir)
    return repo_dir, planned_blocks(repo_dir, "--configfile", "config/full.yaml")


def test_full_config_collect_command_stays_small_and_specs_are_a_file(full_blocks):
    _, blocks = full_blocks
    assert len(blocks["collect_tables"]) == 1
    [command] = [
        line for line in blocks["collect_tables"][0].splitlines() if "--specs" in line
    ]
    assert len(command.encode()) < 100_000
    assert "--specs results/run_specs.json" in command
    assert len(blocks["run_specs"]) == 1


def test_full_config_plans_sim10k_panel_and_nothing_at_sim10k_bfull(full_blocks):
    _, blocks = full_blocks
    assert len(blocks["sim10k_panel"]) == 1
    assert not [
        b for rule_blocks in blocks.values() for b in rule_blocks if "sim10k/Bfull" in b
    ]


def test_bfull_stage_plans_one_pplacer_and_no_epa_ng_dipper_r_or_oracle(full_blocks):
    repo_dir, _ = full_blocks
    blocks = planned_blocks(
        repo_dir,
        "--configfile", "config/full.yaml", "config/full-bfull.yaml",
    )  # fmt: skip
    assert len(blocks["pplacer_place"]) == 1
    for rule in ("epa_ng_place", "dipper_r_place", "oracle_region"):
        assert rule not in blocks
