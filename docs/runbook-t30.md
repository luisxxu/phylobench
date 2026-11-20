# T30 runbook: Bfull\* confirmation on `gg2-wol`

One stage, 37 cells, on `highmem` only. Same conventions as `docs/runbook-t29.md`: run from `barnacle-login`, in `/ddn_scratch/lxxu/phylobench`, `phylobench` env active, inside `tmux`, SLURM account `knightlab`, archive directory `/ddn_scratch/lxxu/phylobench-archive`. One command per block; after each step, paste back what the "Paste" line says. Every `snakemake` run goes through the profile `profiles/slurm-highmem`, targets come before `--configfile`, and the user runs everything; agents never touch the cluster. Never run two `snakemake` invocations in the same working directory. Hand-typed `srun` always carries `-n 1` (see "One task" below). Cluster facts: `docs/cluster.md`.

What T30 asks: is any core challenger as accurate as production SEPP at Bfull\* on the T31 `gg2-wol` set, and what does each arm cost on one CPU model? The queries are the same ones T28 and T29 used at B10k\*, so this is confirmatory in the backbone, not in the queries. Plan and decisions: `TICKETS.md` T30, `design.md` sections 3 and 9. The only placeholder in this file is `<T30_MERGE_SHA>`.

## Pre-registered rules

- **Primary**, per region, core challengers only (`raw/sepp-pinned`, `hmmalign/bscampp`): paired difference in the fraction of queries within 3 edges against `raw/sepp-prod`; cluster bootstrap over held-out groups (a clade is one cluster, a singleton its own), B = 10,000, seed 0; an unplaced query or a DNF challenger is a miss; a DNF reference makes the region undecidable. **Non-inferior** iff the one-sided 95% lower bound is above -0.02.
- **Headline** "as accurate as production SEPP" iff non-inferior at full-length **and** V4-150. Expected joint power about 65% at N = 3,000, stated up front. A challenger with either region `not-shown` is "not shown"; otherwise `undecidable`.
- **Secondary, descriptive:** mean edge-distance difference with its CI, with the 0.25-edge gate printed as undecidable; the G-only skeleton exact-edge difference (k = 0), labelled "descriptive (skeleton exact-edge)", convention in `design.md` section 3; per-bin summary.
- **COST rows** (`hmmalign/{pplacer-scampp, epa-ng-scampp, bscampp-b5000, apples2, dipper-m}` at full-length and V4-150) carry the role `COST (descriptive)`: no verdict, no headline. Their accuracy is never a verdict.
- **Wording:** truth is the G-leaf genome position in GG2. That uDance kept the WoL2 topology fixed is UNVERIFIED. Always carry the apples2 caveat: GG2 truth was seeded by uDance, which used APPLES-2.
- **No upp cell.** T29 S0b was unsized (B2: timeout at 8,402 s for 300 queries, Xeon 6242, `b-cn-25`, peak 31,102 MB) and S2 did not run. The S0b row is printed as text in the cost table; no upp-minus-hmmalign difference exists.
- **Memory stop:** T28 `sepp-prod` full-length peak 96 GiB x 2.9 = 278 GiB, under 1,400 GB. Satisfied.
- **Cost stop: 80 node-hours** (`sacct`). Estimate 25 to 35, UNVERIFIED. A SEPP full-length timeout is a DNF row, not a rerun.

## Matrix

| Item | Value |
|---|---|
| Config | `config/t30.yaml` |
| Profile | `profiles/slurm-highmem` |
| Dataset, backbone | `gg2-wol`, `Bfull` only |
| Regions | T28's frozen 9: `full-length`, `V4-150`, `V3-V4`, `V1-V2`, `V1-V3`, `V4-V5`, `V5-V6`, `V6-V8`, `V7-V9c` |
| Core (27) | `raw/sepp-prod` x 9 (reference), `raw/sepp-pinned` x 9, `hmmalign/bscampp` x 9 |
| COST (10) | `hmmalign/{pplacer-scampp, epa-ng-scampp, bscampp-b5000, apples2, dipper-m}` x {`full-length`, `V4-150`} |
| Not at Bfull\* | `dipper-m-gpu` (T22 GPU row stands), `pplacer` (T22 capped row stands), `epa-ng`, `dipper-r` (B10k\* only) |
| Threads | 16 per job, SLURM and tool flag |

Nothing is reused or rerun: no cell of this config has an earlier run_id (checked by `tests/test_t30_config.py`).

Profile: `hmmalign` and `guppy_tog` (shared) and the 11 place rules (exclusive) run on `highmem` (`b-cn-25`, `b-cn-26`, Intel Xeon Gold 6242), each with `tasks: 1`. `hmmbuild`, `upp_align` and `refprep_me_tree` stay on `long` (not planned here). At most two exclusive jobs run at once, so the stage takes about 1 to 3 days of wall time with queueing.

## One task

On an `--exclusive` 64-CPU node `srun -c 16` without `-n 1` starts 4 tasks (T29 S0b incident). The profile sets `tasks: 1`; the SLURM plugin 2.8.0 turns that into `--ntasks=1` and jobstep 0.6.1 runs `srun -n1`. The offline dry run proves the profile config only. The guard is the `sacct` NTasks check below, run mid-run and after: any step with NTasks other than 1 is a STOP.

## Expected jobs in the dry run

Cluster repo (the T31 outputs are present): **134 jobs**.

| Rule | Jobs |
|---|---|
| `hmmalign`, `mask_hmmalign`, `similarity` | 9 each |
| `sepp_prod`, `sepp_pinned_place`, `bscampp_place` | 9 each |
| `pplacer_scampp_place`, `epa_ng_scampp_place`, `bscampp_b5000_place`, `apples2_place`, `dipper_m_place` | 2 each |
| `dipper_m_input`, `dipper_to_jplace` | 2 each |
| `guppy_tog` | 26 (9 + 9 + 4 x 2) |
| `evaluate` | 37 |
| `run_specs`, `collect_tables`, `all` | 1 each |

Absent (a STOP if any appears): `extract_regions`, `gg2_unzip`, `gg2_wol_build`, `gg2_build`, `backbone_stockholm`, `hmmbuild`, `refprep_*`, `upp_*`, `oracle_region`, `pplacer_place`, `epa_ng_place`, `dipper_r_*`, `dipper_m_gpu_place`, `dipper_gpu_sif`, `concordance`, and anything for B10k or `gg2`.

Greps on the dry-run log: `slurm_partition=highmem` 72, `tasks=1` 72, `--exclusive` 37.

## Cost estimate (UNVERIFIED)

From T22's Bfull\* rows, scaled to this N (full-length 3,000 queries; the 8 amplicon regions 18,937 in all). `sepp-prod` and `sepp-pinned` full-length: about 3 to 3.5 node-hours each (T22: about 650 s fixed + about 3.6 s per query x 3,000); their 8 amplicon cells about 8 to 8.6 node-hours each; `bscampp` x 9 about 2.5; COST x 10 about 4.1; `hmmalign` x 9 and `guppy_tog` x 26 about 1.8; short rules about 3; the skeleton step about 1. Total **about 25 to 35 node-hours** (T22 S2 actual: 51 node-hours for 89 cells). SEPP peaks near 280 GB. Pre-registered stop: **80 node-hours**.

## Preflight

```bash
tmux new -s t30
```

Paste: nothing.

```bash
cd /ddn_scratch/lxxu/phylobench && conda activate phylobench
```

Paste: nothing.

```bash
git fetch origin && git status --short
```

Paste: the output (expect empty).

```bash
git checkout --detach <T30_MERGE_SHA>
```

Paste: the last line of output.

```bash
git rev-parse HEAD
```

Paste: the SHA.

Nothing in the pipeline code may have changed since the T31 prep, or finished prep jobs rerun. Expect only `bscampp.smk` and `upp.smk` listed (T29); anything in `align`, `refprep` or `envs/` is a STOP.

```bash
git diff 49e7565 HEAD --stat -- envs/ workflow/rules/
```

Paste: the output.

```bash
python -c "import phylobench, snakemake; print(phylobench.__file__, snakemake.__version__)"; conda list -n phylobench "snakemake-executor-plugin-slurm*" | grep slurm
```

Paste: the output (expect a path under `/ddn_scratch/lxxu/phylobench/src`, `9.27.0`, plugin `2.8.0`, jobstep `0.6.1`).

```bash
python -m pytest -q tests/test_skeleton.py tests/test_skeleton_table.py tests/test_decision_table.py tests/test_t30_config.py tests/test_profile.py tests/test_summaries.py tests/test_within_k.py tests/test_tables.py
```

Paste: the last line.

The upstream files must exist: 9 region query files, plus the Bfull tree, reference tree and backbone alignment, the held-out table, the full query files, the me-tree, the epa-ng model and the backbone HMM (18 in all).

```bash
for r in full-length V4-150 V3-V4 V1-V2 V1-V3 V4-V5 V5-V6 V6-V8 V7-V9c; do ls -d results/gg2-wol/regions/$r/queries.fasta; done 2>&1 | grep -c "^results"; ls -d data/gg2-wol/Bfull/tree.nwk data/gg2-wol/Bfull/reference_tree.nwk data/gg2-wol/Bfull/aligned-dna-sequences.fasta data/gg2-wol/heldout.tsv data/gg2-wol/queries.fasta data/gg2-wol/queries.aligned.fasta results/ref-prep/gg2-wol/Bfull/me-tree.nwk results/ref-prep/gg2-wol/Bfull/epa-ng.model results/gg2-wol/Bfull/align/backbone.hmm 2>&1 | grep -c "^\(data\|results\)"
```

Paste: the two numbers (expect 9 and 9; otherwise paste the `ls` errors).

## highmem check

```bash
sinfo -p highmem -o "%P %N %c %m %t"
```

Paste: the output (expect `b-cn-25` and `b-cn-26`, 64 CPUs, 1546528 MB).

```bash
for n in b-cn-25 b-cn-26; do srun -A knightlab -p highmem -w $n -n 1 -c 1 --mem=1G -t 00:05:00 lscpu | grep "Model name"; done
```

Paste: the output (expect Intel Xeon Gold 6242 twice).

```bash
squeue -p highmem -h | wc -l
```

Paste: the number (jobs queued or running ahead of us).

```bash
sbatch --test-only -A knightlab -p highmem --exclusive --ntasks=1 --mem=1500000 -t 4320 --wrap "true"
```

Paste: the output (expect a start-time estimate, no error).

## One-task evidence from T29

T29 S1 and the S0 probes ran as one task each. Count steps with NTasks other than 1 (expect 0 bad).

```bash
sacct -P -u lxxu -A knightlab -S 2026-11-18 --format=JobID,NTasks,AllocCPUS,NodeList | awk -F'|' 'NR==1{next} $1 ~ /\.[0-9]+$/ && $2!="" {n++; if ($2!=1) bad++; if ($2!=1) print} END{print n+0, "steps;", bad+0, "NTasks != 1"}'
```

Paste: the output (the void S0b attempt, step `2316728.0`, is expected as the one line with NTasks 4).

## t29-final checksum guard

`collect_tables` rewrites `results/*.tsv` with the planned cells only, so the T29 tables must be archived first. `cp -n` is not used; a mismatch is copied under a new name and recorded.

```bash
mkdir -p logs; test -s /ddn_scratch/lxxu/phylobench-archive/t29-final/accuracy.tsv && test -s /ddn_scratch/lxxu/phylobench-archive/t29-final/runs.tsv && echo "t29-final present" || echo "STOP: t29-final missing from the archive; do not continue"
```

Paste: the output. On `STOP`, do not continue.

```bash
sha256sum results/runs.tsv results/accuracy.tsv results/concordance.tsv; cat /ddn_scratch/lxxu/phylobench-archive/t29-final/SHA256SUMS
```

Paste: the output. If the three `results/` checksums equal the `runs.tsv`, `accuracy.tsv` and `concordance.tsv` lines of `SHA256SUMS`, no copy is needed. If any differs, copy under a new name and paste the checksums:

```bash
for t in runs accuracy concordance; do cp results/$t.tsv /ddn_scratch/lxxu/phylobench-archive/$t-pre-t30.tsv; done; sha256sum /ddn_scratch/lxxu/phylobench-archive/*-pre-t30.tsv
```

Paste: the output (or "no copy needed").

## Dry run

```bash
snakemake all -n --profile profiles/slurm-highmem --configfile config/t30.yaml > logs/t30-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A40 "^Job stats" logs/t30-dryrun.txt | head -50
```

Paste: the output. STOP unless the total is 134 and every rule count matches "Expected jobs in the dry run" above; STOP if any absent rule appears.

```bash
grep -c "slurm_partition=highmem" logs/t30-dryrun.txt; grep -c "tasks=1\(,\|$\)" logs/t30-dryrun.txt; grep -c -- "--exclusive" logs/t30-dryrun.txt
```

Paste: the three numbers (expect 72, 72, 37).

```bash
snakemake all --list-changes code --profile profiles/slurm-highmem --configfile config/t30.yaml; echo "--- params"; snakemake all --list-changes params --profile profiles/slurm-highmem --configfile config/t30.yaml; echo "--- input"; snakemake all --list-changes input --profile profiles/slurm-highmem --configfile config/t30.yaml
```

Paste: the output. Expect `code` empty, `params` only `run_specs`, `input` only `collect_tables`; anything else is a STOP (finished prep would rerun).

Wait for my confirmation before the run.

## Run

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t30.start
```

Paste: nothing.

```bash
snakemake all --profile profiles/slurm-highmem --jobs 6 --keep-going --configfile config/t30.yaml 2>&1 | tee logs/t30.log
```

Paste: nothing until it ends; then the last 30 lines of `logs/t30.log`.

From a second window, check once a day while it runs. Number of `srun` steps and any with NTasks other than 1 (job rows have empty NTasks, so `-X` is not used; `-S` is required):

```bash
sacct -P -r highmem -u lxxu -A knightlab -S "$(cat logs/t30.start)" --format=JobID,NTasks,AllocCPUS,NodeList | awk -F'|' 'NR==1{print; next} $1 ~ /\.[0-9]+$/ {n++; if ($2!=1) bad++} END{print n+0, "srun steps;", bad+0, "NTasks != 1"}'
```

Paste: the last line. **STOP and `scancel` the whole stage if any NTasks is not 1.**

Node-hours so far (conservative: counts every job of yours in the window):

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t30.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number. **STOP and `scancel` the stage above 80.**

## After the run

```bash
grep -c "Error in rule" logs/t30.log
```

Paste: the number, and the matching lines if it is not 0. Two-failure rule: if two jobs of one class `(rule, backbone)` fail, stop the stage and paste; do not rerun. A SEPP timeout or OOM is a DNF row (data), not a failure.

Stage gate.

```bash
python workflow/scripts/check_stage.py config/t30.yaml
```

Paste: the full output. `FAIL` lines block the stage; `INFO` lines are OOM or timeout DNFs (data).

The same NTasks command over the finished stage (expect at most 72 `srun` steps, 0 bad, nodes `b-cn-25` and `b-cn-26`; a job that died before `srun` has no step):

```bash
sacct -P -r highmem -u lxxu -A knightlab -S "$(cat logs/t30.start)" --format=JobID,NTasks,AllocCPUS,NodeList | awk -F'|' 'NR==1{print; next} $1 ~ /\.[0-9]+$/ {n++; if ($2!=1) bad++} END{print n+0, "srun steps;", bad+0, "NTasks != 1"}'
```

Paste: the full output. **STOP if any NTasks is not 1 or a step ran off `b-cn-25` and `b-cn-26`.**

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t30.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number (stop is 80).

`cpu_model` check (`check_stage.py` does not look at it). The `hmmbuild` row has no region and ran on `long`, so it is excluded. Expect 72 rows: 9 `hmmalign` + 37 place-class + 26 `convert`.

```bash
python - <<'EOF'
import pandas as pd
runs = pd.read_csv("results/runs.tsv", sep="\t")
phases = ["align", "place", "convert", "align+place", "align+place+convert"]
rows = runs[(runs["dataset"] == "gg2-wol") & (runs["backbone"] == "Bfull") & runs["region"].notna() & runs["phase"].isin(phases)]
print(len(rows), "rows (expect 72)")
print(rows.groupby(["cpu_model", "threads"], dropna=False).size())
bad = rows[~rows["cpu_model"].fillna("").str.contains("Xeon.*Gold 6242") | (rows["threads"] != 16)]
print("STOP: rows off Xeon Gold 6242 or threads != 16" if len(bad) else "cpu_model ok")
print(bad[["run_id", "phase", "cpu_model", "threads"]].to_string(index=False))
EOF
```

Paste: the output. **STOP if it prints `STOP`.**

## Tables

Primary table: every challenger cell against `raw/sepp-prod` (28 rows: 37 cells less the 9 references).

```bash
python workflow/scripts/within_k_table.py config/t30.yaml | tee results/t30-within3.txt
```

Paste: the full output.

The skeleton step reads the 200k-leaf backbone and every jplace; it runs on a short node, one task.

```bash
srun -A knightlab -p short -n 1 -c 2 --mem=64G -t 04:00:00 python workflow/scripts/skeleton_table.py config/t30.yaml
```

Paste: the output (expect `N scored, M skipped (no accuracy rows)`, N + M = 37; a skipped cell is a DNF).

```bash
python workflow/scripts/within_k_table.py config/t30.yaml --accuracy results/skeleton.tsv -k 0 | tee results/t30-skeleton-within0.txt
```

Paste: the full output. Descriptive (skeleton exact-edge), never a verdict.

```bash
python workflow/scripts/decision_table.py config/t30.yaml | tee results/t30-decision.txt
```

Paste: the full output (the decision table with `COST (descriptive)` rows, then the headline per core arm).

Cost table: end-to-end per cell, the `hmmalign` rows, the one-time `hmmbuild` row (excluded from e2e; printed with its `cpu_model`), and the S0b text line.

```bash
python - <<'EOF' | tee results/t30-cost.txt
import pandas as pd
from phylobench.summaries import end_to_end
runs = pd.read_csv("results/runs.tsv", sep="\t")
runs = runs[(runs["dataset"] == "gg2-wol") & (runs["backbone"] == "Bfull")]
cols = ["region", "aln", "arm", "query_n", "place_s", "e2e_s", "peak_mem_mb", "cpu_model", "e2e_status"]
print(end_to_end(runs)[cols].sort_values(["region", "arm", "aln"]).to_string(index=False))
align = runs[runs["arm"].isna() & runs["region"].notna() & (runs["phase"] == "align")]
print(align[["region", "aln", "wall_s", "peak_mem_mb", "cpu_model", "status"]].sort_values("region").to_string(index=False))
print(runs[runs["arm"].isna() & runs["region"].isna() & (runs["phase"] == "align")][["tool", "wall_s", "peak_mem_mb", "cpu_model", "status"]].to_string(index=False))
print("S0b (T29, upp/bscampp at Bfull*, 300 queries, 16 threads, Xeon Gold 6242, b-cn-25): timeout at 8,402 s, peak 31,102 MB; UNSIZED, no upp cell")
EOF
```

Paste: the full output.

```bash
python - <<'EOF' | tee results/t30-bins.txt
import pandas as pd
from phylobench.summaries import bin_summary
acc = pd.read_csv("results/accuracy.tsv", sep="\t")
acc = acc[acc["run_id"].str.startswith("gg2-wol/Bfull/")]
print(bin_summary(acc).to_string(index=False))
EOF
```

Paste: the output (per-bin results are descriptive; never pool bins across backbone sizes).

## Archive and copy to the laptop

```bash
git rev-parse HEAD > results/GIT_SHA
```

Paste: nothing.

```bash
mkdir -p /ddn_scratch/lxxu/phylobench-archive/t30-final
```

Paste: nothing.

```bash
cp results/runs.tsv results/accuracy.tsv results/concordance.tsv results/skeleton.tsv results/t30-*.txt results/GIT_SHA /ddn_scratch/lxxu/phylobench-archive/t30-final/
```

Paste: nothing.

```bash
tar czf /ddn_scratch/lxxu/phylobench-archive/t30-final/logs.tgz config/t30.yaml profiles/slurm-highmem/config.yaml logs/t30*; echo $?
```

Paste: the number.

```bash
cd /ddn_scratch/lxxu/phylobench-archive/t30-final && sha256sum runs.tsv accuracy.tsv concordance.tsv skeleton.tsv t30-*.txt logs.tgz > SHA256SUMS && cat SHA256SUMS
```

Paste: the output.

On the laptop, from the repo root. The first command must print nothing.

```bash
test -e results/t30-final && echo "STOP: results/t30-final already exists"
```

Paste: the output.

```bash
scp -r lxxu@barnacle.ucsd.edu:/ddn_scratch/lxxu/phylobench-archive/t30-final results/ && cd results/t30-final && shasum -a 256 -c SHA256SUMS
```

Paste: the output.

## Result

Filled in by the close-out PR (`Closes #61`) after the run. Template:

- Run: date, SHA, jobs (of 134), `Error in rule` count, `check_stage.py` result, node-hours (stop 80), wall time, NTasks line (n steps, 0 bad), `cpu_model` check (72 rows), archive `t30-final` checksums verified on the laptop.
- Confirmatory caveat: the queries are the same ones T28 and T29 used at B10k\*; Bfull\* is a new backbone, not a new query set.
- Headline per core arm (`results/t30-decision.txt`): `raw/sepp-pinned` and `hmmalign/bscampp`, full-length and V4-150 within-3 difference and lower bound, and the words "as accurate as production SEPP", "not shown" or "undecidable".
- Per-region table (28 rows): core verdicts; COST rows descriptive only.
- Secondary: mean edge-distance difference with CI and the 0.25 gate (undecidable); skeleton exact-edge difference (descriptive).
- Cost (`results/t30-cost.txt`, Xeon Gold 6242 only): e2e and peak memory at full-length and V4-150 per arm, the `hmmbuild` row, and the S0b line as the upp Bfull\* evidence ("unsized", no upp-minus-hmmalign difference).
- apples2 caveat: GG2 truth was seeded by uDance, which used APPLES-2. Truth is the G-leaf genome position in GG2; uDance WoL2-topology fixedness is UNVERIFIED.
- DNF rows, if any, listed as data.
