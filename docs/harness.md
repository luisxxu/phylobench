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
- On Linux `ru_maxrss` (`maxrss_mb`) is inherited across fork/exec, so it is at least the harness's own RSS at spawn (small for the standalone CLI). Never call `run()` in-process from a large parent when `maxrss_mb` matters.

## 5. Status caveats

- `oom` means exit by SIGKILL (returncode -9 or 137). An external SIGKILL not caused by the OOM killer is mislabeled `oom` on the psutil path; `returncode` is kept for inspection.
- The SLURM wall limit must exceed `--timeout` plus margin, otherwise SLURM kills the harness before it writes JSON (T10 uses 86400).

## 6. GPU memory

`--gpu` sums `used_memory` over all compute processes on the visible GPUs (`nvidia-smi --query-compute-apps`). That is correct for an exclusive `--gres` allocation and works with `apptainer --nv`.

## 7. Arm contract (T8)

Adding a tool arm means adding files only: `workflow/rules/<tool>.smk`, `config/arms/<arm>.yaml`, `envs/<tool>.yaml`, `docs/tools/<tool>.md`. T8's `Snakefile`, `evaluate.smk`, `smoke.smk`, `config/schema.yaml`, `config/smoke.yaml` and `ci.yml` are never edited. The template is `workflow/rules/pplacer.smk` and `config/arms/pplacer.yaml`. `*-gpu` arms are manual tier, excluded from the smoke matrix.

**Arm YAML** (strict, validated by `config/schema.yaml`; the file stem is the arm name):

| Key | Meaning |
|---|---|
| `tool` | required; the `tool` column of the arm's phases in `runs.tsv` (the `convert` phase is always `guppy`) |
| `aligned` | required; `true` runs over every `alns` entry in the config, `false` runs once with `aln=raw` on the unaligned region queries |
| `phases` | required; ordered map phase name to its `run.json` filename in the cell directory. The first entry is the gating phase that produced `placements.jplace`. Listing `convert: placements.tog.run.json` triggers the shared `guppy_tog` rule |
| `smoke_min_exact` | optional; minimum exact placements on the full-length cell for each aln in the smoke test |

**Cell** = `results/{dataset}/{backbone}/{region}/{aln}/{arm}/`; `run_id` is `dataset/backbone/region/aln/arm`. The rule writes `placements.jplace` (jplace v3 on the backbone tree, edge numbers as in `data/{dataset}/{backbone}/tree.nwk`) and one `run.json` per listed phase. `evaluate` adds `accuracy.tsv` and `evaluate.run.json`. EMP cells (`emp-10k`, `emp-100k`; region `V4-150`, aln `hmmalign` only) have no truth tree, so the `concordance` rule adds `concordance.tsv` and `concordance.run.json` instead (the `evaluate` phase row in `runs.tsv`); see [`docs/datasets/emp.md`](datasets/emp.md).

**Inputs** a rule may use: aligned queries `results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta`; raw queries `results/{dataset}/regions/{region}/queries.fasta`; backbone `data/{dataset}/{backbone}/{tree.nwk,aligned-dna-sequences.fasta,raxml-info.txt}`; truth `data/{dataset}/{backbone}/reference_tree.nwk` (the reference tree with the queries attached; T8 defines this contract and T18 must produce it); ref-prep outputs under `results/ref-prep/{dataset}/{backbone}/`. For EMP datasets the backbone files, `me-tree.nwk`, `epa-ng.model` and `backbone.hmm` are links to `data/gg2`, `results/ref-prep/gg2` and `results/gg2` (no `reference_tree.nwk`), so the gg2 stage must run before EMP.

**Run contract**: every tool call is `{MEASURE} --out <run.json> --threads {threads} --tool-version '<pin>' --env-lock envs/<env>.linux-64.pin.txt -- <cmd>` followed by the DNF touch snippet of section 2; `threads` is `config["threads"]` and is also passed to the tool's own flag; `conda:` is the tool env; `log:` is set.

**DNF mapping**: in `runs.tsv` a harness status other than `ok` becomes `status=DNF` with `dnf_reason` set to the harness status (`oom|timeout|error`). `evaluate` reads only the gating `run.json`: not ok gives a header-only `accuracy.tsv` and the jplace is never opened; ok with an empty jplace is a hard failure. A failing `evaluate` itself fails the job (a scoring bug is never a DNF row).

**Concordance** (EMP): `results/concordance.tsv` has one row per query per EMP cell, `edge_dist` against the `sepp-prod` placement on the same backbone, `label=concordance`, NaN where either side left the query unplaced (never skipna). The file is header-only for a cell when either the arm's or `sepp-prod`'s gating run is not ok; ok with an empty jplace is a hard failure. The `sepp-prod` cells score against themselves (all 0); consumers filter those self rows out (`run_id == reference_run_id`). `results/accuracy.tsv` never contains concordance rows.

**Forced DNF**: `smoke.smk` runs `sleep 30` under `--timeout 1` into the cell `.../forced-dnf/` (`forced_dnf` in `config/smoke.yaml`), yielding a `place` row with `status=DNF`, `dnf_reason=timeout` and no accuracy rows.

**Acceptance checks** (`workflow/scripts/check_smoke.py`, CI per arm, logic in `phylobench.tables.smoke_failures`): every aligned-query file has the backbone's width; every planned `(run_id, phase)` appears exactly once and nothing else; the forced row is DNF/timeout; every other row is `ok` with non-null `version` and `env_lock_sha256`; every non-forced cell has one non-NaN accuracy row per fixture query; arms with `smoke_min_exact` meet it on the full-length cell. CI also fails if Snakemake logs `Failed to install conda environment from pin file` (silent fallback to the unlocked YAML).
