import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import psutil
import pytest

from phylobench import measure

SRC = str(Path(__file__).resolve().parent.parent / "src")
CLI_ENV = {**os.environ, "PYTHONPATH": SRC}
RUN_KEYS = {
    "status",
    "returncode",
    "timeout_s",
    "wall_s",
    "cpu_s",
    "peak_mem_mb",
    "mem_source",
    "psutil_peak_mb",
    "cgroup_peak_mb",
    "maxrss_mb",
    "node",
    "cpu_model",
    "threads",
    "version",
    "env_lock_sha256",
    "command",
    "start_utc",
}


def py(code):
    return [sys.executable, "-c", code]


def gone(pids, within=5):
    procs = []
    for pid in pids:
        try:
            procs.append(psutil.Process(pid))
        except psutil.NoSuchProcess:
            pass
    _, alive = psutil.wait_procs(procs, timeout=within)
    return [p for p in alive if p.status() != psutil.STATUS_ZOMBIE] == []


def read_pids(path):
    return [int(x) for x in Path(path).read_text().split()]


@pytest.mark.parametrize(
    ("timed_out", "returncode", "expected"),
    [
        (True, 0, "timeout"),
        (True, -9, "timeout"),
        (False, 0, "ok"),
        (False, -9, "oom"),
        (False, 137, "oom"),
        (False, 1, "error"),
        (False, -11, "error"),
        (False, -15, "error"),
    ],
)
def test_classify_status(timed_out, returncode, expected):
    assert measure.classify_status(timed_out, returncode) == expected


def test_parse_cgroup_path():
    assert (
        measure.parse_cgroup_path("0::/user.slice/a.scope\n") == "/user.slice/a.scope"
    )
    assert measure.parse_cgroup_path("12:memory:/foo\n11:cpu:/bar\n") is None


def test_parse_gpu_used_mb():
    assert measure.parse_gpu_used_mb("123\n456\n[N/A]\n") == 579
    assert measure.parse_gpu_used_mb("") == 0


def test_maxrss_mib():
    assert measure.maxrss_mib(1048576, "darwin") == 1.0
    assert measure.maxrss_mib(1024, "linux") == 1.0


def test_select_mem_source():
    assert measure.select_mem_source(None) == "psutil"
    assert measure.select_mem_source(12.0) == "cgroup"


def test_read_cgroup_peak(tmp_path):
    assert measure.read_cgroup_peak(tmp_path) is None
    (tmp_path / "memory.peak").write_text("4096\n")
    assert measure.read_cgroup_peak(tmp_path) == 4096
    (tmp_path / "memory.peak").write_text("max\n")
    assert measure.read_cgroup_peak(tmp_path) is None


def test_run_ok_schema():
    result = measure.run(
        ["sh", "-c", "exit 0"], threads=4, version="1.2", env_lock_sha256="ab"
    )
    assert set(result) == RUN_KEYS
    assert result["status"] == "ok"
    assert result["returncode"] == 0
    assert result["mem_source"] in ("cgroup", "psutil")
    assert result["threads"] == 4
    assert result["version"] == "1.2"
    assert result["command"] == ["sh", "-c", "exit 0"]
    for key in ("wall_s", "cpu_s", "psutil_peak_mb", "maxrss_mb", "peak_mem_mb"):
        assert isinstance(result[key], float)
    assert isinstance(result["node"], str)
    assert result["start_utc"].endswith("+00:00")


def test_run_gpu_keys(monkeypatch):
    monkeypatch.setattr(measure, "sample_gpu_used_mb", lambda: 7.0)
    monkeypatch.setattr(measure, "read_gpu_model", lambda: "FakeGPU")
    result = measure.run(["sh", "-c", "exit 0"], gpu=True)
    assert result["gpu_peak_mb"] == 7.0
    assert result["gpu_model"] == "FakeGPU"


def test_run_exit_code_is_error():
    result = measure.run(py("import sys; sys.exit(3)"))
    assert (result["status"], result["returncode"]) == ("error", 3)


def test_run_self_sigkill_is_oom():
    result = measure.run(py("import os; os.kill(os.getpid(), 9)"))
    assert (result["status"], result["returncode"]) == ("oom", -9)


def test_run_exit_137_is_oom():
    assert measure.run(["sh", "-c", "exit 137"])["status"] == "oom"


def test_run_timeout_kills_whole_tree(tmp_path):
    pids = tmp_path / "pids"
    script = f"sleep 60 & echo $! >> {pids}; sleep 60 & echo $! >> {pids}; wait"
    result = measure.run(["sh", "-c", script], timeout=1)
    assert result["status"] == "timeout"
    assert result["wall_s"] < 5
    assert gone(read_pids(pids))


def test_run_sweeps_orphans(tmp_path):
    pid_file = tmp_path / "pid"
    result = measure.run(["sh", "-c", f"sleep 60 & echo $! > {pid_file}"])
    assert result["status"] == "ok"
    assert gone(read_pids(pid_file))


FORKING_PARENT = """
import os, subprocess, sys, time
ready = sys.argv[1:]
child = (
    "import os,sys,time; b=os.urandom(50*2**20);"
    " open(sys.argv[1],'w').close(); time.sleep(2.5)"
)
procs = [subprocess.Popen([sys.executable, "-c", child, r]) for r in ready]
while not all(os.path.exists(r) for r in ready):
    time.sleep(0.05)
time.sleep(1)
for p in procs:
    p.wait()
"""


def test_run_forking_peak_exceeds_single_process(tmp_path):
    script = tmp_path / "parent.py"
    script.write_text(FORKING_PARENT)
    ready = [str(tmp_path / f"ready{i}") for i in range(3)]
    out = tmp_path / "run.json"
    # Fresh process: Linux ru_maxrss is inherited across fork/exec, so run()
    # called from a large pytest process would report the parent's RSS.
    done = cli("--out", str(out), "--", sys.executable, str(script), *ready)
    assert done.returncode == 0
    result = json.loads(out.read_text())
    assert result["status"] == "ok"
    assert result["psutil_peak_mb"] >= 150
    assert result["psutil_peak_mb"] > result["maxrss_mb"]


def test_run_cpu_time():
    busy = "import time; t=time.time(); exec('while time.time()-t<1: pass')"
    result = measure.run(["sh", "-c", '"$0" -c "$1"', sys.executable, busy])
    assert result["cpu_s"] >= 0.5


def cli(*args, **kwargs):
    return subprocess.run(
        [sys.executable, "-m", "phylobench.measure", *args],
        env=kwargs.pop("env", CLI_ENV),
        capture_output=True,
        text=True,
        **kwargs,
    )


def test_cli_writes_json(tmp_path):
    out = tmp_path / "run.json"
    done = cli("--out", str(out), "--threads", "16", "--", "sh", "-c", "echo hi")
    assert done.returncode == 0
    assert done.stdout == "hi\n"
    data = json.loads(out.read_text())
    assert data["status"] == "ok"
    assert data["threads"] == 16
    assert not (tmp_path / "run.json.tmp").exists()


def test_cli_env_lock_hash(tmp_path):
    lock = tmp_path / "lock.txt"
    lock.write_text("pkg=1.0\n")
    out = tmp_path / "run.json"
    done = cli(
        "--out", str(out), "--env-lock", str(lock), "--tool-version", "9", "--", "true"
    )
    assert done.returncode == 0
    data = json.loads(out.read_text())
    assert data["env_lock_sha256"] == hashlib.sha256(lock.read_bytes()).hexdigest()
    assert data["version"] == "9"


def test_cli_nonzero_tool_still_exit_zero(tmp_path):
    out = tmp_path / "run.json"
    done = cli("--out", str(out), "--", "sh", "-c", "exit 5")
    assert done.returncode == 0
    assert json.loads(out.read_text())["status"] == "error"


@pytest.mark.parametrize(
    "extra",
    [
        ["--"],
        [],
        ["--env-lock", "/nonexistent/lock", "--", "true"],
        ["--", "/nonexistent/binary"],
    ],
)
def test_cli_config_errors_exit_2(tmp_path, extra):
    out = tmp_path / "run.json"
    assert cli("--out", str(out), *extra).returncode == 2
    assert not out.exists()


def test_cli_missing_out_dir_exits_2(tmp_path):
    out = tmp_path / "nodir" / "run.json"
    assert cli("--out", str(out), "--", "true").returncode == 2


def test_cli_gpu_without_nvidia_smi_exits_2(tmp_path):
    out = tmp_path / "run.json"
    env = {**CLI_ENV, "PATH": ""}
    assert cli("--out", str(out), "--gpu", "--", "true", env=env).returncode == 2
    assert not out.exists()


def test_cli_sigterm_kills_child_no_json(tmp_path):
    out = tmp_path / "run.json"
    pid_file = tmp_path / "pid"
    proc = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "phylobench.measure",
            "--out",
            str(out),
            "--",
            "sh",
            "-c",
            f"echo $$ > {pid_file}; sleep 60",
        ],
        env=CLI_ENV,
    )
    deadline = time.monotonic() + 10
    while not pid_file.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    time.sleep(0.3)
    proc.send_signal(signal.SIGTERM)
    assert proc.wait(timeout=10) == 143
    assert gone(read_pids(pid_file))
    assert not out.exists()


def test_console_script(tmp_path):
    if shutil.which("phylobench-measure") is None:
        if os.environ.get("CI"):
            pytest.fail("phylobench-measure not on PATH after pip install -e .")
        pytest.skip("not installed: run pip install -e . --no-deps")
    out = tmp_path / "run.json"
    done = subprocess.run(["phylobench-measure", "--out", str(out), "--", "true"])
    assert done.returncode == 0
    assert json.loads(out.read_text())["status"] == "ok"


def test_real_cgroup_peak():
    cgroup_dir = measure.own_cgroup_dir()
    if cgroup_dir is None or measure.read_cgroup_peak(cgroup_dir) is None:
        pytest.skip("cgroup v2 memory.peak not readable")
    result = measure.run(py("import os,time; b=os.urandom(50*2**20); time.sleep(1.5)"))
    assert result["psutil_peak_mb"] > 0
    cgroup_peak = result["cgroup_peak_mb"]
    assert cgroup_peak is None or cgroup_peak >= result["psutil_peak_mb"] - 10
