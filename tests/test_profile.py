"""Dry-run the SLURM profile with the local executor; nothing is submitted."""

import re

import yaml
from emp_helpers import ROOT
from test_workflow_dag import planned_blocks, repo  # noqa: F401

PROFILE = ROOT / "profiles" / "slurm" / "config.yaml"
THREADED_TOOL_RULES = {
    "hmmbuild",
    "hmmalign",
    "refprep_me_tree",
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
    "guppy_tog",
}
GPU_RULE = "dipper_m_gpu_place"
NODE_MEM_MB = 514000


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
    assert threaded == THREADED_TOOL_RULES | {GPU_RULE}
    for rule in THREADED_TOOL_RULES:
        for job in blocks[rule]:
            res = resources(job)
            assert res["slurm_partition"] == "long", rule
            assert res["runtime"] == "4320", rule
            assert "--exclusive" in res["slurm_extra"], rule
            assert int(res["mem_mb"]) <= NODE_MEM_MB, rule
    for job in blocks[GPU_RULE]:
        res = resources(job)
        assert res["slurm_partition"] == "gpu" and res["gpu"] == "1"
        assert "--nodelist=b-gpu-21" in res["slurm_extra"]
        assert "--exclusive" in res["slurm_extra"]
        assert "gres" not in res
        assert int(res["mem_mb"]) <= 508000
    for rule in ("evaluate", "collect_tables"):
        for job in blocks[rule]:
            res = resources(job)
            assert res["slurm_partition"] == "short"
            assert "--exclusive" not in res["slurm_extra"]
    for jobs in blocks.values():
        assert all("--no-requeue" in resources(job)["slurm_extra"] for job in jobs)
