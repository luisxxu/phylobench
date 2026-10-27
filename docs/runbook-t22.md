# T22 runbook: full benchmark on Barnacle

Stages S1, S2, S3 and the archive. Run everything from `barnacle-login`, in `/ddn_scratch/lxxu/phylobench`, with the `phylobench` env active, inside `tmux`. SLURM account `knightlab`. Archive directory `/ddn_scratch/lxxu/phylobench-archive`. One command per block; after each step, paste back what the "Paste" line says.

Never run two `snakemake` invocations in the same working directory (including the P0 probe: wait for it to finish). Cluster facts: `docs/cluster.md`. Measurement contract: `docs/harness.md`.

## Matrix

`config/full.yaml` is the whole matrix (304 cells); the overlays pick a backbone. What `full.yaml` excludes, and why, is in its header.

| Stage | Command configs | Profile | Cells |
|---|---|---|---|
| S1 B10k\* | `full.yaml` + `full-b10k.yaml` | `profiles/slurm-pilot` | 215 (gg2 171, sim10k 22, emp-10k 11, emp-100k 11); the 33 pilot cells are already done, so 182 are new |
| S2 Bfull\* | `full.yaml` + `full-bfull.yaml` | `profiles/slurm` | 89 (gg2 73 including the one capped pplacer cell, emp-10k 8, emp-100k 8) |
| S3 final collect | `full.yaml` alone | `profiles/slurm` | 0 new jobs |

New place-rule jobs per stage (the dry run may show fewer where a cell is already complete, e.g. the pilot or the P0 probe cell; it must never show more):

| Rule | S1 (new) | S2 |
|---|---|---|
| `sepp_prod` | 10 | 11 |
| `sepp_pinned_place` | 10 | 11 |
| `pplacer_place` | 19 | 1 |
| `pplacer_scampp_place` | 19 | 11 |
| `epa_ng_scampp_place` | 19 | 11 |
| `bscampp_place` | 19 | 11 |
| `epa_ng_place` | 19 | 0 |
| `apples2_place` | 19 | 11 |
| `dipper_m_place` | 19 | 11 |
| `dipper_m_gpu_place` | 19 | 11 |
| `dipper_r_place` | 10 | 0 |

S2 sends 22 jobs to `highmem` (`sepp_prod`, `sepp_pinned_place`, 11 each); the other place rules run on whole `long` nodes. S1's B10k\* peaks set this split: SEPP full-length reaches 103 GB and scaled ~2.9x at Bfull\* in P0, too close to a 514 GB node; the SCAMPP-family arms peak at 49 GB or less. S1 uses `profiles/slurm-pilot`, which keeps every rule on `short`.

Projected cost (pilot walls scaled): S1 about 17 node-hours (ceiling 30), S2 about 365 (ceiling 550). If S2's `sacct` total passes its ceiling, stop and paste.

`dipper-m` at Bfull\* stays planned in `config/full.yaml` until the P0 probe says otherwise. To drop it, add `- {backbone: [Bfull], arm: [dipper-m, dipper-m-gpu]}` to `exclude` in `config/full.yaml` (the commented line at its end) and update the counts in `tests/test_full_config.py`; the stage tables above then lose the two dipper-m rows' S2 jobs.

## Preflight

Start a persistent session.

```bash
tmux new -s t22
```

Paste: nothing.

```bash
cd /ddn_scratch/lxxu/phylobench
```

Paste: nothing.

```bash
conda activate phylobench
```

Paste: nothing.

```bash
git fetch origin
```

Paste: nothing.

```bash
git status --short
```

Paste: the output (expect empty).

Check out the T22 merge commit, detached, so the run is tied to one SHA.

```bash
git checkout --detach <T22_MERGE_SHA>
```

Paste: the last line of output.

```bash
git rev-parse HEAD
```

Paste: the SHA.

`envs/phylobench.yaml` must be unchanged since the pilot, or the pilot's finished jobs rerun.

```bash
git diff 7fde884 HEAD --stat -- envs/phylobench.yaml
```

Paste: the output (expect empty).

```bash
python -c "import phylobench, snakemake; print(phylobench.__file__, snakemake.__version__)"
```

Paste: the output (expect a path under `/ddn_scratch/lxxu/phylobench/src` and `9.27.0`).

```bash
python -m pytest -q tests/test_full_config.py tests/test_gates.py
```

Paste: the last line.

```bash
mkdir -p /ddn_scratch/lxxu/phylobench-archive logs
```

Paste: nothing.

Highmem and queue checks.

```bash
sinfo -p highmem -N -o "%N %T %m %c"
```

Paste: the output.

```bash
scontrol show partition highmem
```

Paste: the output.

```bash
scontrol show node b-cn-25
```

Paste: the output.

```bash
squeue -p highmem
```

Paste: the output.

```bash
sbatch --test-only -A knightlab -p highmem --exclusive --mem=1500000 -t 4320 --wrap true
```

Paste: the output.

```bash
sbatch --test-only -A knightlab -p long --exclusive --mem=500000 -t 4320 --wrap true
```

Paste: the output.

```bash
sbatch --test-only -A knightlab -p short --exclusive --mem=500000 -t 1440 --wrap true
```

Paste: the output.

CPU model of the node classes (recorded per row in `runs.tsv`, but read it once here).

```bash
srun -A knightlab -p highmem -n 1 -t 5 lscpu | grep "Model name"
```

Paste: the output.

```bash
srun -A knightlab -p long -n 1 -t 5 lscpu | grep "Model name"
```

Paste: the output.

## S1: B10k\* (215 cells)

Dry run to a file.

```bash
snakemake all -n --profile profiles/slurm-pilot --configfile config/full.yaml config/full-b10k.yaml > logs/t22-s1-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A60 "^Job stats" logs/t22-s1-dryrun.txt | head -70
```

Paste: the output. Expect `0` for `gg2_download`, `gg2_unzip`, `gg2_build`, `emp_download`, `dipper_gpu_sif` (anything else there means data on scratch is missing or changed), one `sim10k_panel`, one `run_specs`, one `collect_tables`, and the place-rule counts of the S1 column at most.

Which finished jobs would rerun because their code changed (expect only `collect_tables`):

```bash
snakemake all --list-changes code --profile profiles/slurm-pilot --configfile config/full.yaml config/full-b10k.yaml
```

Paste: the output.

Wait for my confirmation of account, partitions and node-hours before submitting. Then record the start time.

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t22-s1.start
```

Paste: nothing.

Submit.

```bash
snakemake all --profile profiles/slurm-pilot --jobs 6 --keep-going --configfile config/full.yaml config/full-b10k.yaml 2>&1 | tee logs/t22-s1.log
```

Paste: nothing until it ends; then the last 30 lines of `logs/t22-s1.log`. If you detach from tmux (`Ctrl-b d`), reattach with `tmux attach -t t22`.

Monitor from a second tmux window (`Ctrl-b c`).

```bash
squeue -u lxxu -o "%.10i %.28j %.8T %.10M %.12P %R"
```

Paste: the output when asked, or when something looks stuck.

```bash
grep -c "Error in rule" logs/t22-s1.log
```

Paste: the number, and the matching lines if it is not 0.

After it ends, run the stage gate (same configs as the run).

```bash
python workflow/scripts/check_stage.py config/full.yaml config/full-b10k.yaml
```

Paste: the full output. `FAIL` lines block the stage. `WARN dnf_error` lines need your acknowledgement. `INFO` lines are OOM or timeout DNFs (data, not failures).

```bash
cp results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-t22-s1.tsv
```

Paste: nothing.

Node-hours since the S1 start (other jobs of yours in that window are counted too).

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t22-s1.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number.

Gate for S2: `check_stage.py` has no `FAIL` and node-hours are under 30.

## S2: Bfull\* (89 cells)

Long-queue check (exclusive `long` and `highmem` jobs wait for whole nodes).

```bash
squeue -p long -t PD -h | wc -l
```

Paste: the number.

```bash
squeue -p highmem -h | wc -l
```

Paste: the number.

Dry run to a file.

```bash
snakemake all -n --profile profiles/slurm --configfile config/full.yaml config/full-bfull.yaml > logs/t22-s2-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A60 "^Job stats" logs/t22-s2-dryrun.txt | head -70
```

Paste: the output. Expect `0` for `gg2_download`, `gg2_unzip`, `gg2_build`, `emp_download`, `dipper_gpu_sif`, `hmmbuild`, `refprep_me_tree`, `backbone_stockholm` (Bfull\* prep is already done), and the S2 column counts at most, with `pplacer_place` exactly 1 and no `epa_ng_place`, `dipper_r_place` or `oracle_region`.

```bash
grep -c "slurm_partition=highmem" logs/t22-s2-dryrun.txt
```

Paste: the number (expect 22, fewer only for finished cells).

```bash
snakemake all --list-changes code --profile profiles/slurm --configfile config/full.yaml config/full-bfull.yaml
```

Paste: the output (expect only `collect_tables`).

Wait for my confirmation, then record the start time.

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t22-s2.start
```

Paste: nothing.

Submit.

```bash
snakemake all --profile profiles/slurm --keep-going --configfile config/full.yaml config/full-bfull.yaml 2>&1 | tee logs/t22-s2.log
```

Paste: nothing until it ends; then the last 30 lines of `logs/t22-s2.log`. Expect several days. A `dipper_m_gpu_place` job pending for more than 24 hours on `b-gpu-21`: paste `squeue` and stop to ask.

```bash
squeue -u lxxu -o "%.10i %.28j %.8T %.10M %.12P %R"
```

Paste: the output when asked.

```bash
grep -c "Error in rule" logs/t22-s2.log
```

Paste: the number, and the matching lines if it is not 0.

After it ends:

```bash
python workflow/scripts/check_stage.py config/full.yaml config/full-bfull.yaml
```

Paste: the full output.

```bash
cp results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-t22-s2.tsv
```

Paste: nothing.

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t22-s2.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number (ceiling 550).

Gate for S3: no `FAIL`, node-hours under the ceiling.

## S3: final collect (`full.yaml` alone)

This plans all 304 cells, so a missing `run.json` from S1 or S2 would submit that job on `profiles/slurm`. The dry run must show only collection jobs before you submit.

```bash
snakemake all -n --profile profiles/slurm --configfile config/full.yaml > logs/t22-s3-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A12 "^Job stats" logs/t22-s3-dryrun.txt
```

Paste: the output. Expect only `all`, `run_specs`, `collect_tables`, and at most one `emp_manifest`. Anything else: stop, do not submit.

Submit.

```bash
snakemake all --profile profiles/slurm --configfile config/full.yaml 2>&1 | tee logs/t22-s3.log
```

Paste: the last 15 lines.

The T22 acceptance check (every planned cell has a row; scored cells have exactly `query_n` rows; gates).

```bash
python workflow/scripts/check_stage.py config/full.yaml
```

Paste: the full output.

## Archive

```bash
git rev-parse HEAD > results/GIT_SHA
```

Paste: nothing.

```bash
du -sh results
```

Paste: the output.

Per-cell working directories, temp dirs and regenerable intermediates are left out; the placements (`placements.jplace`), `run.json` files, per-cell tables and logs stay.

```bash
tar czf /ddn_scratch/lxxu/phylobench-archive/t22-final.tgz --exclude='work' --exclude='tmp0' --exclude='*.in.fasta' --exclude='*.tog.tre' --exclude='*.afa' --exclude='*.aligned.fasta' results config profiles envs data/gg2/MANIFEST.tsv data/emp/MANIFEST.tsv logs
```

Paste: nothing yet; the next step shows its exit status.

```bash
echo $?
```

Paste: the number.

```bash
cd /ddn_scratch/lxxu/phylobench-archive
```

Paste: nothing.

```bash
sha256sum t22-final.tgz > SHA256SUMS
```

Paste: nothing.

```bash
mkdir t22-final
```

Paste: nothing.

```bash
cp /ddn_scratch/lxxu/phylobench/results/runs.tsv /ddn_scratch/lxxu/phylobench/results/accuracy.tsv /ddn_scratch/lxxu/phylobench/results/concordance.tsv t22-final/
```

Paste: nothing.

```bash
cd t22-final
```

Paste: nothing.

```bash
sha256sum runs.tsv accuracy.tsv concordance.tsv > SHA256SUMS
```

Paste: nothing.

```bash
cat SHA256SUMS
```

Paste: the output.

```bash
ls -l /ddn_scratch/lxxu/phylobench-archive
```

Paste: the output.

Copy the three tables to the laptop (`concordance.tsv` is about 150 MB), from the repo root there, and verify.

```bash
scp -r barnacle-login:/ddn_scratch/lxxu/phylobench-archive/t22-final results/
```

Paste: nothing.

```bash
cd results/t22-final
```

Paste: nothing.

```bash
shasum -a 256 -c SHA256SUMS
```

Paste: the output.

## Two-failure rule

A job class is `(rule, backbone)`, e.g. `(bscampp_place, Bfull)`. A failure is Snakemake's `Error in rule` (a SLURM kill, a missing `run.json`, a scoring bug); a harness DNF (`oom`, `timeout`, `error` in `runs.tsv`) is data, not a failure.

- First failure of a class: paste the `Error in rule` block and its log; I diagnose and fix, then you resubmit the same stage command with `--rerun-incomplete` added.
- Second failure of the same class: stop the stage and paste. No third attempt.
- If the tmux session or the submitting Snakemake dies: do not restart; paste `squeue -u lxxu` first (running SLURM jobs keep running, and a restart would double-submit them).
