"""Dry-run the SLURM profile with the local executor; nothing is submitted."""

import re

import yaml
from emp_helpers import ROOT
from test_workflow_dag import planned_blocks, repo  # noqa: F401

PROFILE = ROOT / "profiles" / "slurm" / "config.yaml"
PREP_RULES = {"hmmbuild", "hmmalign", "refprep_me_tree", "guppy_tog"}
PLACE_RULES = {
    "sepp_prod",
    "sepp_pinned_place",
    "pplacer_place",
    "pplacer_scampp_place",
    "epa_ng_scampp_place",
    "bscampp_place",
    "epa_ng_place",
    "apples2_place",
    "dipper_m_place",
    "dipper_r_place",
}
GPU_RULE = "dipper_m_gpu_place"
NODE_MEM_MB = 514000
GPU_MAX_MEM_MB = 16 * 12288  # gpu partition MaxMemPerCPU x 16 cpus/task


def resources(block: str) -> dict[str, str]:
    line = re.search(r"resources: (.*)", block)[1]
    return dict(re.findall(r"(\w+)=(.*?)(?:, (?=\w+=)|$)", line))


def test_profile_never_retries_and_never_groups():
    profile = yaml.safe_load(PROFILE.read_text())
    assert profile["retries"] == 0
    assert not {k for k in profile if "group" in k}


def test_slurm_resources_by_rule(repo):  # noqa: F811
    blocks = planned_blocks(
        repo, "-c", "2", "--profile", "profiles/slurm", "--executor", "local",
        "--configfile", "config/smoke.yaml",
    )  # fmt: skip
    threaded = {
        rule
        for rule, jobs in blocks.items()
        if any("--threads 2" in job for job in jobs)
    }
    assert threaded == PREP_RULES | PLACE_RULES | {GPU_RULE}
    for rule in PREP_RULES | PLACE_RULES:
        for job in blocks[rule]:
            res = resources(job)
            assert res["slurm_partition"] == "long", rule
            assert res["runtime"] == "4320", rule
            assert int(res["mem_mb"]) <= NODE_MEM_MB, rule
            assert ("--exclusive" in res["slurm_extra"]) == (rule in PLACE_RULES), rule
    for job in blocks[GPU_RULE]:
        res = resources(job)
        assert res["slurm_partition"] == "gpu" and res["gpu"] == "1"
        assert "--nodelist=b-gpu-21" in res["slurm_extra"]
        assert "--exclusive" in res["slurm_extra"]
        assert "gres" not in res
        assert int(res["mem_mb"]) <= GPU_MAX_MEM_MB
    for rule in ("evaluate", "collect_tables"):
        for job in blocks[rule]:
            res = resources(job)
            assert res["slurm_partition"] == "short"
            assert "--exclusive" not in res["slurm_extra"]
    for jobs in blocks.values():
        assert all("--no-requeue" in resources(job)["slurm_extra"] for job in jobs)


def test_pilot_profile_only_moves_long_rules_to_short():
    profile = yaml.safe_load(PROFILE.read_text())
    pilot = yaml.safe_load((ROOT / "profiles/slurm-pilot/config.yaml").read_text())
    for rule in PREP_RULES | PLACE_RULES:
        res = pilot["set-resources"][rule]
        assert (res["slurm_partition"], res["runtime"]) == ("short", 1440), rule
        profile["set-resources"][rule] = res
    assert pilot == profile
