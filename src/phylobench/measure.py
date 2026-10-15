import argparse
import contextlib
import hashlib
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import threading
import time
from datetime import UTC, datetime
from pathlib import Path

import psutil

SAMPLE_INTERVAL_S = 0.5
CGROUP_ROOT = Path("/sys/fs/cgroup")
MIB = 1024 * 1024


def classify_status(timed_out: bool, returncode: int) -> str:
    if timed_out:
        return "timeout"
    if returncode == 0:
        return "ok"
    if returncode in (-signal.SIGKILL, 128 + signal.SIGKILL):
        return "oom"
    return "error"


def parse_cgroup_path(text: str) -> str | None:
    for line in text.splitlines():
        if line.startswith("0::"):
            return line[3:]
    return None


def parse_gpu_used_mb(text: str) -> float:
    total = 0.0
    for line in text.splitlines():
        with contextlib.suppress(ValueError):
            total += float(line.strip())
    return total


def maxrss_mib(ru_maxrss: int, platform_name: str) -> float:
    # ru_maxrss is KiB on Linux but bytes on macOS
    divisor = MIB if platform_name == "darwin" else 1024
    return ru_maxrss / divisor


def select_mem_source(cgroup_peak_mb: float | None) -> str:
    return "psutil" if cgroup_peak_mb is None else "cgroup"


def read_cgroup_peak(cgroup_dir: Path) -> int | None:
    try:
        return int((cgroup_dir / "memory.peak").read_text())
    except (OSError, ValueError):
        return None


def own_cgroup_dir() -> Path | None:
    try:
        path = parse_cgroup_path(Path("/proc/self/cgroup").read_text())
    except OSError:
        return None
    return None if path is None else CGROUP_ROOT / path.lstrip("/")


def sample_tree_rss(root: psutil.Process) -> int:
    total = 0
    try:
        procs = [root, *root.children(recursive=True)]
    except psutil.Error:
        return 0
    for p in procs:
        with contextlib.suppress(psutil.Error):
            total += p.memory_info().rss
    return total


def sample_gpu_used_mb() -> float | None:
    try:
        out = subprocess.run(
            [
                "nvidia-smi",
                "--query-compute-apps=used_memory",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=5,
            check=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return parse_gpu_used_mb(out.stdout)


def read_cpu_model() -> str | None:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return None


def read_gpu_model() -> str | None:
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    names = out.stdout.strip().splitlines()
    return names[0].strip() if names else None


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def kill_tree(root: psutil.Process, pgid: int) -> None:
    try:
        descendants = root.children(recursive=True)
    except psutil.Error:
        descendants = []
    with contextlib.suppress(ProcessLookupError, PermissionError):
        os.killpg(pgid, signal.SIGKILL)
    for p in descendants:
        with contextlib.suppress(psutil.Error):
            p.kill()
    with contextlib.suppress(psutil.Error):
        root.kill()


def write_json(path: Path, data: dict) -> None:
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    os.replace(tmp, path)


def run(
    cmd: list[str],
    *,
    timeout: float | None = None,
    gpu: bool = False,
    threads: int | None = None,
    version: str | None = None,
    env_lock_sha256: str | None = None,
) -> dict:
    cgroup_dir = own_cgroup_dir()
    peak_start = read_cgroup_peak(cgroup_dir) if cgroup_dir else None
    start_utc = datetime.now(UTC).isoformat()
    t0 = time.monotonic()
    proc = subprocess.Popen(cmd, start_new_session=True)
    root = psutil.Process(proc.pid)
    deadline = None if timeout is None else t0 + timeout

    stop = threading.Event()
    state = {"rss": 0, "gpu": 0.0, "timed_out": False}

    def sampler():
        while True:
            state["rss"] = max(state["rss"], sample_tree_rss(root))
            if gpu:
                used = sample_gpu_used_mb()
                if used is not None:
                    state["gpu"] = max(state["gpu"], used)
            if deadline is not None and time.monotonic() > deadline:
                state["timed_out"] = True
                kill_tree(root, proc.pid)
                return
            if stop.wait(SAMPLE_INTERVAL_S):
                return

    thread = threading.Thread(target=sampler, daemon=True)
    thread.start()
    try:
        _, status, rusage = os.wait4(proc.pid, 0)
    except BaseException:
        kill_tree(root, proc.pid)
        with contextlib.suppress(ChildProcessError):
            os.wait4(proc.pid, 0)
        raise
    finally:
        stop.set()
        thread.join()
        with contextlib.suppress(ProcessLookupError, PermissionError):
            os.killpg(proc.pid, signal.SIGKILL)
    wall_s = time.monotonic() - t0
    returncode = os.waitstatus_to_exitcode(status)
    proc.returncode = returncode

    peak_end = read_cgroup_peak(cgroup_dir) if cgroup_dir else None
    cgroup_peak_mb = (
        peak_end / MIB
        if peak_start is not None and peak_end is not None and peak_end > peak_start
        else None
    )
    psutil_peak_mb = state["rss"] / MIB
    mem_source = select_mem_source(cgroup_peak_mb)
    result = {
        "status": classify_status(state["timed_out"], returncode),
        "returncode": returncode,
        "timeout_s": timeout,
        "wall_s": wall_s,
        "cpu_s": rusage.ru_utime + rusage.ru_stime,
        "peak_mem_mb": cgroup_peak_mb if mem_source == "cgroup" else psutil_peak_mb,
        "mem_source": mem_source,
        "psutil_peak_mb": psutil_peak_mb,
        "cgroup_peak_mb": cgroup_peak_mb,
        "maxrss_mb": maxrss_mib(rusage.ru_maxrss, sys.platform),
        "node": platform.node(),
        "cpu_model": read_cpu_model(),
        "threads": threads,
        "version": version,
        "env_lock_sha256": env_lock_sha256,
        "command": list(cmd),
        "start_utc": start_utc,
    }
    if gpu:
        result["gpu_peak_mb"] = state["gpu"]
        result["gpu_model"] = read_gpu_model()
    return result


def _exit_on_sigterm(signum, frame):
    raise SystemExit(128 + signum)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="phylobench-measure")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--tool-version")
    parser.add_argument("--env-lock", type=Path)
    parser.add_argument("cmd", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    cmd = args.cmd[1:] if args.cmd[:1] == ["--"] else args.cmd

    def config_error(message: str) -> int:
        print(f"phylobench-measure: {message}", file=sys.stderr)
        return 2

    if not cmd:
        return config_error("no command given after --")
    if not args.out.parent.is_dir():
        return config_error(f"--out directory missing: {args.out.parent}")
    if args.gpu and shutil.which("nvidia-smi") is None:
        return config_error("--gpu requires nvidia-smi on PATH")
    env_lock_sha256 = None
    if args.env_lock is not None:
        try:
            env_lock_sha256 = sha256_file(args.env_lock)
        except OSError as e:
            return config_error(f"--env-lock unreadable: {e}")

    signal.signal(signal.SIGTERM, _exit_on_sigterm)
    try:
        result = run(
            cmd,
            timeout=args.timeout,
            gpu=args.gpu,
            threads=args.threads,
            version=args.tool_version,
            env_lock_sha256=env_lock_sha256,
        )
    except (FileNotFoundError, PermissionError) as e:
        return config_error(f"cannot execute {cmd[0]}: {e}")
    except KeyboardInterrupt:
        return 130
    write_json(args.out, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
