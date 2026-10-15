# Measurement harness (`phylobench-measure`)

Every tool call in the workflow is wrapped by `phylobench-measure` (T4). It runs one command, passes stdout/stderr through, and writes a flat `run.json`.

```sh
phylobench-measure --out run.json [--timeout S] [--gpu] [--threads N] [--tool-version STR] [--env-lock PATH] -- CMD...
```

Exit codes: `0` whenever `run.json` was written (whatever the tool's status); `2` pre-start configuration error, no JSON (empty CMD, `--out` parent missing, unreadable `--env-lock`, `--gpu` without `nvidia-smi`, CMD not executable); `130`/`143` on SIGINT/SIGTERM to the harness (process tree killed, no JSON); `1` internal error.

`run.json` fields (memory in MiB, 2^20): `status` (`ok|oom|timeout|error`), `returncode` (negative = signal), `timeout_s`, `wall_s`, `cpu_s` (`os.wait4` user+sys), `peak_mem_mb`, `mem_source` (`cgroup|psutil`), `psutil_peak_mb`, `cgroup_peak_mb` (null unless the own cgroup's `memory.peak` rose during the run), `maxrss_mb` (`os.wait4` `ru_maxrss`), `node`, `cpu_model`, `threads`, `version`, `env_lock_sha256`, `command`, `start_utc`; with `--gpu` also `gpu_peak_mb` and `gpu_model`.

## 1. Invocation from rules

Rules call `f"{sys.executable} -m phylobench.measure"`, resolved once at Snakefile parse time. This is the absolute path of the phylobench env's Python, so it works under `--use-conda`; CMD itself is resolved on the tool env's PATH. For GPU tools on the host, CMD is `apptainer exec --nv ...`. T8 defines the harness invocation once as a Snakefile constant.

Job-side Snakemake must run from the same env (UNVERIFIED, decided in T21).

If the package or psutil is missing the job fails with `ModuleNotFoundError`, exit 1, no JSON. That is a loud job failure, never a DNF.

## 2. DNF contract

The harness exits 0 whenever it wrote `run.json`, so the tool's failure does not fail the job by itself. Only when `run.json` exists and `status != ok` does the rule shell touch the declared tool outputs the tool failed to produce, so the job succeeds and `run.json` survives:

```sh
{harness} --out {output.run_json} -- {tool_cmd}
if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
    for f in {output.tool_files}; do [ -e "$f" ] || : > "$f"; done
fi
```

Consumers read `run.json` first: `status != ok` gives a `status=DNF` row and the outputs are ignored; `status == ok` with a missing or empty output is a hard failure. The T8 smoke includes one forced-DNF case (`--timeout 1`).

## 3. Provenance assertions

T8 asserts that every non-forced row has `status == ok` and non-null `version` and `env_lock_sha256`.

## 4. Memory caveats

- `memory.peak` of the harness's own cgroup includes the harness, the shell, and Snakemake ancestors (about 100-200 MiB) and page cache. It is meaningful only with one rule per SLURM job. T27/T21 decide which source is primary.
- The psutil RSS sum double-counts shared pages, and 0.5 s sampling misses short spikes.
- Both values are always recorded; `peak_mem_mb` follows `mem_source`.

## 5. Status caveats

- `oom` means exit by SIGKILL (returncode -9 or 137). An external SIGKILL not caused by the OOM killer is mislabeled `oom` on the psutil path; `returncode` is kept for inspection.
- The SLURM wall limit must exceed `--timeout` plus margin, otherwise SLURM kills the harness before it writes JSON (T10 uses 86400).

## 6. GPU memory

`--gpu` sums `used_memory` over all compute processes on the visible GPUs (`nvidia-smi --query-compute-apps`). That is correct for an exclusive `--gres` allocation and works with `apptainer --nv`.
