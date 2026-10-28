# T28 runbook: column-window regions on Barnacle

Stages S0 (exclusion-count check), S1 (B10k\*), S2 (Bfull\*), collect. Same conventions as `docs/runbook-t22.md`: run from `barnacle-login`, in `/ddn_scratch/lxxu/phylobench`, `phylobench` env active, inside `tmux`, SLURM account `knightlab`, archive directory `/ddn_scratch/lxxu/phylobench-archive`. One command per block; after each step, paste back what the "Paste" line says. Never run two `snakemake` invocations in the same working directory. Cluster facts: `docs/cluster.md`.

What the column-window regions are: `docs/gg2.md` ("Column-window regions (T28)"). Pre-registered window rule, tolerance and reporting criteria: `TICKETS.md` T28.

## Matrix

`config/t28.yaml` plans only new cells: gg2, the 8 `c` regions (`V4-150c`, `V3-V4c`, `V1-V2c`, `V1-V3c`, `V4-V5c`, `V5-V6c`, `V6-V8c`, `V7-V9c`), arms `sepp-prod`, `sepp-pinned`, `bscampp` (hmmalign), both backbones. The T22 rows (full-length and the primer-matched regions) are reused, nothing there reruns.

| Stage | Command configs | Profile | Cells |
|---|---|---|---|
| S0 check | none (extract only) | none | 0 |
| S1 B10k\* | `t28.yaml` + `t28-b10k.yaml` | `profiles/slurm-pilot` | 24 (8 per arm) |
| S2 Bfull\* | `t28.yaml` + `t28-bfull.yaml` | `profiles/slurm` | 24 (8 per arm) |

Expected jobs per stage (a dry run may show fewer, never more): `extract_regions` 0 (S0 already cut the 8 `c` regions; they are dataset-level, shared by both stages), `hmmalign` 8, `mask_hmmalign` 8, `similarity` 8, `sepp_prod` 8, `sepp_pinned_place` 8, `bscampp_place` 8, `guppy_tog` 16, `evaluate` 24. Never `oracle_region`, `pplacer_place`, `epa_ng_place`, `dipper_*`, `apples2_place`.

Cost estimate (unverified, scaled from the T22 rows: 8 regions instead of 11 cells per arm): S1 about 1 to 2 node-hours (ceiling 5), S2 about 12 to 16 (ceiling 25). If S2's `sacct` total passes its ceiling, stop and paste.

`collect_tables` rewrites `results/runs.tsv`, `accuracy.tsv`, `concordance.tsv` with the planned cells of the configs only, so the T22 tables are copied aside first (preflight).

## Preflight

```bash
tmux new -s t28
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

```bash
git checkout --detach <T28_MERGE_SHA>
```

Paste: the last line of output.

```bash
git rev-parse HEAD
```

Paste: the SHA.

```bash
git diff 7fde884 HEAD --stat -- envs/phylobench.yaml
```

Paste: the output (expect empty; otherwise finished jobs rerun).

```bash
python -c "import phylobench, snakemake; print(phylobench.__file__, snakemake.__version__)"
```

Paste: the output (expect a path under `/ddn_scratch/lxxu/phylobench/src` and `9.27.0`).

```bash
python -m pytest -q tests/test_column_regions.py tests/test_t28_config.py tests/test_full_config.py
```

Paste: the last line.

Keep the T22 tables. The check block must print `T22 tables present`; if it prints `STOP`, do not continue and paste the output (a missing T22 table means the wrong directory or a lost result).

```bash
mkdir -p /ddn_scratch/lxxu/phylobench-archive logs
```

Paste: nothing.

```bash
for t in runs accuracy concordance; do test -s results/$t.tsv || echo "STOP: results/$t.tsv missing or empty"; done; test -s results/runs.tsv -a -s results/accuracy.tsv -a -s results/concordance.tsv && echo "T22 tables present: $(wc -l < results/runs.tsv) runs.tsv lines, $(grep -c -E '(V4-150|V3-V4|V1-V2|V1-V3|V4-V5|V5-V6|V6-V8|V7-V9)c' results/runs.tsv) lines with a c region"
```

Paste: the output (expect `T22 tables present` and 0 lines with a c region; a nonzero count means a T28 stage already rewrote the tables, so do not archive them as T22).

`cp -n` never overwrites an existing archive file, so re-running this block after a stage or a tmux crash cannot replace the T22 copy with a T28 table.

```bash
cp -n results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-t22-final.tsv; cp -n results/accuracy.tsv /ddn_scratch/lxxu/phylobench-archive/accuracy-t22-final.tsv; cp -n results/concordance.tsv /ddn_scratch/lxxu/phylobench-archive/concordance-t22-final.tsv
```

Paste: nothing.

```bash
ls -l /ddn_scratch/lxxu/phylobench-archive/*-t22-final.tsv; sha256sum /ddn_scratch/lxxu/phylobench-archive/*-t22-final.tsv
```

Paste: the output. If a T22 archive file already existed before this step, it is the T22 copy: keep it.

## S0: exclusion counts (no placement jobs)

First the primer-matched regions from T22 (already on scratch). The UMI explanation is unverified until these numbers are read.

```bash
for r in V4-150 V3-V4 V1-V2 V1-V3 V4-V5 V5-V6 V6-V8 V7-V9; do echo "== $r"; cat results/gg2/regions/$r/exclusion_counts.tsv; done
```

Paste: the output.

How many excluded sources are MJ0* UMI amplicons, per region.

```bash
for r in V1-V2 V1-V3 V7-V9 V3-V4; do echo "$r $(cut -f1 results/gg2/regions/$r/exclusions.tsv | grep -c '^MJ0') of $(wc -l < results/gg2/regions/$r/exclusions.tsv) excluded"; done
```

Paste: the output.

Now cut the 8 column-window regions (8 short local jobs, seconds each; reads only `data/gg2/queries*.fasta`).

```bash
snakemake -c 2 --configfile config/t28.yaml results/gg2/regions/V4-150c/queries.fasta results/gg2/regions/V3-V4c/queries.fasta results/gg2/regions/V1-V2c/queries.fasta results/gg2/regions/V1-V3c/queries.fasta results/gg2/regions/V4-V5c/queries.fasta results/gg2/regions/V5-V6c/queries.fasta results/gg2/regions/V6-V8c/queries.fasta results/gg2/regions/V7-V9c/queries.fasta 2>&1 | tail -5
```

Paste: the output.

```bash
for r in V4-150c V3-V4c V1-V2c V1-V3c V4-V5c V5-V6c V6-V8c V7-V9c; do echo "== $r"; cat results/gg2/regions/$r/extract_regions.log; cat results/gg2/regions/$r/exclusion_counts.tsv; done
```

Paste: the output (window per region and the counts). Expect `extracted` near 1,000 for every `c` region and the original-region counts above for comparison. A region with `extracted` under 100 stops the plan: paste and wait.

## S1: B10k\* (24 cells)

Dry run to a file.

```bash
snakemake all -n --profile profiles/slurm-pilot --configfile config/t28.yaml config/t28-b10k.yaml > logs/t28-s1-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A30 "^Job stats" logs/t28-s1-dryrun.txt | head -40
```

Paste: the output. Expect `0` (absent) for `gg2_*`, `hmmbuild`, `refprep_*`, `backbone_stockholm`, `oracle_region`, and the counts of the Matrix section at most.

```bash
snakemake all --list-changes code --profile profiles/slurm-pilot --configfile config/t28.yaml config/t28-b10k.yaml
```

Paste: the output (expect only `collect_tables`).

Wait for my confirmation, then record the start time.

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t28-s1.start
```

Paste: nothing.

```bash
snakemake all --profile profiles/slurm-pilot --jobs 6 --keep-going --configfile config/t28.yaml config/t28-b10k.yaml 2>&1 | tee logs/t28-s1.log
```

Paste: nothing until it ends; then the last 30 lines of `logs/t28-s1.log`.

```bash
grep -c "Error in rule" logs/t28-s1.log
```

Paste: the number, and the matching lines if it is not 0.

Stage gate.

```bash
python workflow/scripts/check_stage.py config/t28.yaml config/t28-b10k.yaml
```

Paste: the full output. `FAIL` lines block the stage; `INFO` lines are OOM or timeout DNFs (data).

```bash
cp results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-t28-s1.tsv
```

Paste: nothing.

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t28-s1.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number. Gate for S2: no `FAIL` and under 5 node-hours.

## S2: Bfull\* (24 cells)

```bash
snakemake all -n --profile profiles/slurm --configfile config/t28.yaml config/t28-bfull.yaml > logs/t28-s2-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A30 "^Job stats" logs/t28-s2-dryrun.txt | head -40
```

Paste: the output (same expectations as S1).

```bash
grep -c "slurm_partition=highmem" logs/t28-s2-dryrun.txt
```

Paste: the number (expect 16: `sepp_prod` and `sepp_pinned_place`, 8 each; `bscampp_place` runs on `long`).

```bash
snakemake all --list-changes code --profile profiles/slurm --configfile config/t28.yaml config/t28-bfull.yaml
```

Paste: the output (expect only `collect_tables`).

Wait for my confirmation, then record the start time.

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t28-s2.start
```

Paste: nothing.

```bash
snakemake all --profile profiles/slurm --keep-going --configfile config/t28.yaml config/t28-bfull.yaml 2>&1 | tee logs/t28-s2.log
```

Paste: nothing until it ends; then the last 30 lines of `logs/t28-s2.log`.

```bash
grep -c "Error in rule" logs/t28-s2.log
```

Paste: the number, and the matching lines if it is not 0.

```bash
python workflow/scripts/check_stage.py config/t28.yaml config/t28-bfull.yaml
```

Paste: the full output.

```bash
cp results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-t28-s2.tsv
```

Paste: nothing.

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t28-s2.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number (ceiling 25).

## Collect (both backbones)

All 48 cells; every job is already done, so the dry run must show only collection jobs.

```bash
snakemake all -n --profile profiles/slurm --configfile config/t28.yaml > logs/t28-collect-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A12 "^Job stats" logs/t28-collect-dryrun.txt
```

Paste: the output (expect only `all`, `run_specs`, `collect_tables`). Anything else: stop, do not submit.

```bash
snakemake all --profile profiles/slurm --configfile config/t28.yaml 2>&1 | tee logs/t28-collect.log
```

Paste: the last 15 lines.

```bash
python workflow/scripts/check_stage.py config/t28.yaml
```

Paste: the full output.

## Archive and copy to the laptop

```bash
git rev-parse HEAD > results/GIT_SHA
```

Paste: nothing.

```bash
mkdir -p /ddn_scratch/lxxu/phylobench-archive/t28-final
```

Paste: nothing.

```bash
cp results/runs.tsv results/accuracy.tsv results/concordance.tsv results/GIT_SHA /ddn_scratch/lxxu/phylobench-archive/t28-final/
```

Paste: nothing.

The exclusion tables and window logs are part of the evidence.

```bash
tar czf /ddn_scratch/lxxu/phylobench-archive/t28-final/column-regions.tgz results/gg2/regions/*c config/t28.yaml config/t28-b10k.yaml config/t28-bfull.yaml logs/t28-*
```

Paste: nothing yet; the next step shows its exit status.

```bash
echo $?
```

Paste: the number.

```bash
cd /ddn_scratch/lxxu/phylobench-archive/t28-final
```

Paste: nothing.

```bash
sha256sum runs.tsv accuracy.tsv concordance.tsv column-regions.tgz > SHA256SUMS
```

Paste: nothing.

```bash
cat SHA256SUMS
```

Paste: the output.

On the laptop, from the repo root. The T28 tables go to their own directory `results/t28-final`; the verified T22 tables and `SHA256SUMS` in `results/` (and `results/t22-final`) are never touched. The first command must print nothing.

```bash
test -e results/t28-final && echo "STOP: results/t28-final already exists"
```

Paste: the output.

```bash
scp -r lxxu@barnacle.ucsd.edu:/ddn_scratch/lxxu/phylobench-archive/t28-final results/
```

Paste: nothing.

```bash
cd results/t28-final
```

Paste: nothing.

```bash
shasum -a 256 -c SHA256SUMS
```

Paste: the output.

## Two-failure rule

Same as `docs/runbook-t22.md`: first failure of a job class `(rule, backbone)` is diagnosed and resubmitted with `--rerun-incomplete`; the second stops the stage; a harness DNF (`oom`, `timeout`, `error`) is data, not a failure. If the tmux session dies, paste `squeue -u lxxu` before restarting anything.
