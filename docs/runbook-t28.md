# T28 runbook: `gg2-wol` at B10k\* on Barnacle

Stages S0 (extract regions, apply the frozen region rule), S1 (B10k\* run), collect. Same conventions as `docs/runbook-t22.md`: run from `barnacle-login`, in `/ddn_scratch/lxxu/phylobench`, `phylobench` env active, inside `tmux`, SLURM account `knightlab`, archive directory `/ddn_scratch/lxxu/phylobench-archive`. One command per block; after each step, paste back what the "Paste" line says. Never run two `snakemake` invocations in the same working directory. Cluster facts: `docs/cluster.md`.

What the set is: `docs/gg2.md` ("WoL genome held-out set (T31)"), already built with ref-prep and `hmmbuild` done for both backbones. Decision rule: `design.md` section 3 and `TICKETS.md` T28 ("T31 amendment"). Column-window regions: `docs/gg2.md` ("Column-window regions (T28)"). T22's `gg2` set is not run; T22 numbers stay as published.

## Matrix

`config/t28.yaml`: dataset `gg2-wol`, backbone `B10k`, 9 regions, 5 arms, 45 cells.

| Item | Value |
|---|---|
| Regions (frozen by S0) | `full-length`, `V4-150`, `V3-V4`, `V1-V2`, `V1-V3`, `V4-V5`, `V5-V6`, `V6-V8`, `V7-V9c` |
| Region rule | `full-length` and `V4-150` always; every other primer region with >= 1,000 queries; else its `c` twin if that has >= 1,000 (primer wins) |
| Arms | `sepp-prod` (raw), `sepp-pinned` (raw), `bscampp`, `pplacer-scampp`, `epa-ng-scampp` (all hmmalign) |
| Profile | `profiles/slurm-pilot` |

The expected list comes from the probe counts (primer-matched queries per region: full-length 3,000, V4-150 2,703, V3-V4 2,542, V1-V2 2,119, V1-V3 2,058, V4-V5 2,622, V5-V6 1,913, V6-V8 2,354, V7-V9 139). Only V7-V9 falls under 1,000, so its `c` twin replaces it; the `V7-V9c` count is UNVERIFIED (probe estimate about 2,600) and S0 measures it. If it is under 1,000, S0 stops and `config/t28.yaml` changes before S1.

Profile: `profiles/slurm-pilot`, as T22's S1 used for B10k\*. Its header explains why: on `long` even shared jobs waited days, `short` outranks it on the same nodes, and B10k\* cells take minutes. Place rules get a whole node (`--exclusive`, 500 GB), enough for SEPP's 103 GB full-length peak at B10k\* (T22). `profiles/slurm` (highmem routing) is for Bfull\* and is not used here.

Expected jobs in the S1 dry run (it may show fewer, never more): `hmmalign` 9, `mask_hmmalign` 9, `similarity` 9, `sepp_prod` 9, `sepp_pinned_place` 9, `bscampp_place` 9, `pplacer_scampp_place` 9, `epa_ng_scampp_place` 9, `guppy_tog` 36, `evaluate` 45, `run_specs` 1, `collect_tables` 1, `all` 1. `extract_regions` is 0 (S0 cut them). Absent: `gg2_*`, `backbone_stockholm`, `hmmbuild`, `refprep_*`, `oracle_region`, `pplacer_place`, `epa_ng_place`, `dipper_*`, `apples2_place`, and anything for `gg2` or `Bfull`.

Cost estimate (UNVERIFIED, scaled from T22 S1's 8.1 node-hours for 215 cells at 1,000 queries): 3 to 8 node-hours. Stop and paste if `sacct` passes 15.

`collect_tables` rewrites `results/runs.tsv`, `accuracy.tsv`, `concordance.tsv` with the planned cells of the config only, so the T22 tables are copied aside first (preflight).

## Preflight

```bash
tmux new -s t28
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
python -m pytest -q tests/test_within_k.py tests/test_within_k_table.py tests/test_t28_config.py tests/test_column_regions.py
```

Paste: the last line.

```bash
ls data/gg2-wol/queries.fasta data/gg2-wol/queries.aligned.fasta data/gg2-wol/heldout.tsv data/gg2-wol/B10k/tree.nwk results/gg2-wol/B10k/align/backbone.hmm results/ref-prep/gg2-wol/B10k/me-tree.nwk
```

Paste: the output (expect all six listed, no `No such file`).

Keep the T22 tables. The check block must print `T22 tables present`; if it prints `STOP`, do not continue and paste the output.

```bash
mkdir -p /ddn_scratch/lxxu/phylobench-archive logs
```

Paste: nothing.

```bash
for t in runs accuracy concordance; do test -s results/$t.tsv || echo "STOP: results/$t.tsv missing or empty"; done; test -s results/runs.tsv -a -s results/accuracy.tsv -a -s results/concordance.tsv && echo "T22 tables present: $(wc -l < results/runs.tsv) runs.tsv lines, $(grep -c 'gg2-wol/' results/runs.tsv) gg2-wol lines"
```

Paste: the output (expect `T22 tables present` and 0 gg2-wol lines; a nonzero count means a T28 stage already rewrote the tables, so do not archive them as T22).

`cp -n` never overwrites an existing archive file, so re-running this block cannot replace the T22 copy with a T28 table.

```bash
cp -n results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-t22-final.tsv; cp -n results/accuracy.tsv /ddn_scratch/lxxu/phylobench-archive/accuracy-t22-final.tsv; cp -n results/concordance.tsv /ddn_scratch/lxxu/phylobench-archive/concordance-t22-final.tsv
```

Paste: nothing.

```bash
ls -l /ddn_scratch/lxxu/phylobench-archive/*-t22-final.tsv; sha256sum /ddn_scratch/lxxu/phylobench-archive/*-t22-final.tsv
```

Paste: the output. If a T22 archive file already existed before this step, it is the T22 copy: keep it.

## S0: extract regions and apply the region rule (no placement jobs)

Cuts the 9 expected regions plus the primer-matched `V7-V9` (for its count) from `data/gg2-wol/queries*.fasta`. Short jobs, seconds each. Targets go before `--configfile` (which takes several files and would swallow them), and the run uses the profile: a run without `--use-conda` records a different software environment, so the next profiled run plans `gg2_wol_build` again (seen on 2026-11-10; the rebuild was byte-identical to the T31 checksums). The `c` extraction reads the shipped aligned rows `data/gg2-wol/queries.aligned.fasta`.

```bash
snakemake results/gg2-wol/regions/full-length/queries.fasta results/gg2-wol/regions/V4-150/queries.fasta results/gg2-wol/regions/V3-V4/queries.fasta results/gg2-wol/regions/V1-V2/queries.fasta results/gg2-wol/regions/V1-V3/queries.fasta results/gg2-wol/regions/V4-V5/queries.fasta results/gg2-wol/regions/V5-V6/queries.fasta results/gg2-wol/regions/V6-V8/queries.fasta results/gg2-wol/regions/V7-V9/queries.fasta results/gg2-wol/regions/V7-V9c/queries.fasta -n --profile profiles/slurm-pilot --configfile config/t28.yaml > logs/t28-s0-dryrun.txt 2>&1; grep -A8 "^Job stats" logs/t28-s0-dryrun.txt
```

Paste: the output. Expect only `extract_regions` (10 jobs, or fewer if some exist). If `gg2_wol_build` or any `refprep_*` appears, stop and paste: the data or its provenance changed since T31.

```bash
snakemake results/gg2-wol/regions/full-length/queries.fasta results/gg2-wol/regions/V4-150/queries.fasta results/gg2-wol/regions/V3-V4/queries.fasta results/gg2-wol/regions/V1-V2/queries.fasta results/gg2-wol/regions/V1-V3/queries.fasta results/gg2-wol/regions/V4-V5/queries.fasta results/gg2-wol/regions/V5-V6/queries.fasta results/gg2-wol/regions/V6-V8/queries.fasta results/gg2-wol/regions/V7-V9/queries.fasta results/gg2-wol/regions/V7-V9c/queries.fasta --profile profiles/slurm-pilot --configfile config/t28.yaml 2>&1 | tail -5
```

Paste: the output.

Exclusion counts per region (the primer-matched counts next to `V7-V9c`, with its window).

```bash
for r in full-length V4-150 V3-V4 V1-V2 V1-V3 V4-V5 V5-V6 V6-V8 V7-V9 V7-V9c; do echo "== $r"; cat results/gg2-wol/regions/$r/extract_regions.log; cat results/gg2-wol/regions/$r/exclusion_counts.tsv; done
```

Paste: the output.

Apply the frozen rule. It prints the per-region query counts, the resulting region list, and `STOP` if the list is not the expected one.

```bash
python - <<'EOF'
import csv
MIN = 1000
primer = ["V3-V4", "V1-V2", "V1-V3", "V4-V5", "V5-V6", "V6-V8", "V7-V9"]
expected = ["full-length", "V4-150", "V3-V4", "V1-V2", "V1-V3", "V4-V5", "V5-V6", "V6-V8", "V7-V9c"]

def extracted(region):
    with open(f"results/gg2-wol/regions/{region}/exclusion_counts.tsv") as f:
        return int(dict(csv.reader(f, delimiter="\t"))["extracted"])

counts = {r: extracted(r) for r in ["full-length", "V4-150", *primer, "V7-V9c"]}
print(counts)
frozen = ["full-length", "V4-150"]
for r in primer:
    frozen.append(r if counts[r] >= MIN else (r + "c" if counts.get(r + "c", 0) >= MIN else "NONE:" + r))
print("frozen list:", frozen)
print("OK, matches config/t28.yaml" if frozen == expected else "STOP: list differs from config/t28.yaml (V7-V9c < 1,000, or another primer region < 1,000); do not run S1")
EOF
```

Paste: the output. On `STOP`, paste and wait: the region list and `config/t28.yaml` change in a new commit first.

Wait for my confirmation, then S1.

## S1: B10k\* (45 cells)

Dry run to a file.

```bash
snakemake all -n --profile profiles/slurm-pilot --configfile config/t28.yaml > logs/t28-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A30 "^Job stats" logs/t28-dryrun.txt | head -40
```

Paste: the output (expectations in the Matrix section).

```bash
snakemake all --list-changes code --profile profiles/slurm-pilot --configfile config/t28.yaml
```

Paste: the output (expect only `collect_tables`).

Wait for my confirmation, then record the start time.

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t28.start
```

Paste: nothing.

```bash
snakemake all --profile profiles/slurm-pilot --jobs 6 --keep-going --configfile config/t28.yaml 2>&1 | tee logs/t28.log
```

Paste: nothing until it ends; then the last 30 lines of `logs/t28.log`.

```bash
grep -c "Error in rule" logs/t28.log
```

Paste: the number, and the matching lines if it is not 0.

Stage gate.

```bash
python workflow/scripts/check_stage.py config/t28.yaml
```

Paste: the full output. `FAIL` lines block the stage; `INFO` lines are OOM or timeout DNFs (data).

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t28.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number (stop and paste above 15).

## Collect and within-3 table

`collect_tables` ran as part of `all`; the tables hold the 45 cells only. The within-3 table is the T28 decision statistic: per region, the paired difference in the fraction of queries within 3 edges of the truth (challenger minus `sepp-prod`), cluster bootstrap over held-out groups (B = 10,000, seed 0), one-sided 95% lower bound, `non-inferior` when the bound is above -0.02; an unplaced query (either arm) is a miss; `undecidable` when the `sepp-prod` cell is DNF (no rows).

```bash
python workflow/scripts/within_k_table.py config/t28.yaml | tee results/t28-within3.txt
```

Paste: the full output (45 rows less the 9 `sepp-prod` reference cells, so 36).

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
cp results/runs.tsv results/accuracy.tsv results/concordance.tsv results/t28-within3.txt results/GIT_SHA /ddn_scratch/lxxu/phylobench-archive/t28-final/
```

Paste: nothing.

The exclusion tables and window logs are part of the evidence.

```bash
tar czf /ddn_scratch/lxxu/phylobench-archive/t28-final/regions.tgz results/gg2-wol/regions config/t28.yaml logs/t28*
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
sha256sum runs.tsv accuracy.tsv concordance.tsv t28-within3.txt regions.tgz > SHA256SUMS && cat SHA256SUMS
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

### Result (2026-11-12, `9ad3bc2`)

S0: every region's extracted plus excluded counts sum to 3,000. Extracted: full-length 3,000, V4-150 2,703, V3-V4 2,542, V1-V2 2,119, V1-V3 2,058, V4-V5 2,622, V5-V6 1,913, V6-V8 2,354, V7-V9 139, `V7-V9c` 2,626 (window columns 1167 to 1552, length 376; 374 below the 0.9 coverage tolerance). Frozen list = `config/t28.yaml` (T30 uses it).

Incident: the first S0 run used `-c 2` without the profile and rebuilt `data/gg2-wol` on the login node (software-environment trigger). All 11 data files matched the T31 sha256s; the stale provenance was cleared with `--cleanup-metadata` on the build and region outputs and the existing ref-prep outputs re-stamped with `--touch`; nothing was recomputed. The S0 block above is the corrected command.

S1: 156 of 156 jobs, 0 `Error in rule`, `check_stage.py` all PASS (no DNFs), 4.5 node-hours (estimate 3 to 8), 10.2 h wall. `sepp-prod` full-length peak 98,257 MB (96 GiB, cgroup, AMD EPYC 7302); x 2.9 = 278 GiB, under the T30 1,400 GB stop.

Within-3 table (`results/t28-within3.txt`; diff = challenger minus `sepp-prod` fraction within 3 edges, lower = one-sided 95% bound, cluster bootstrap B = 10,000; NI = non-inferior, lower > -0.02):

| Region (n, clusters) | bscampp | epa-ng-scampp | pplacer-scampp | sepp-pinned |
|---|---|---|---|---|
| full-length (3,000, 2,258) | -0.030 (-0.052) not shown | -0.023 (-0.044) not shown | -0.023 (-0.044) not shown | 0.000 (0.000) NI |
| V4-150 (2,703, 2,084) | +0.004 (-0.011) NI | +0.005 (-0.010) NI | +0.005 (-0.009) NI | +0.000 (-0.001) NI |
| V3-V4 (2,542, 1,972) | -0.001 (-0.015) NI | -0.001 (-0.015) NI | -0.003 (-0.016) NI | 0.000 (0.000) NI |
| V1-V2 (2,119, 1,690) | -0.004 (-0.022) not shown | +0.004 (-0.014) NI | +0.003 (-0.013) NI | 0.000 (0.000) NI |
| V1-V3 (2,058, 1,639) | +0.003 (-0.013) NI | +0.005 (-0.011) NI | +0.008 (-0.007) NI | 0.000 (0.000) NI |
| V4-V5 (2,622, 2,017) | +0.001 (-0.017) NI | -0.000 (-0.018) NI | -0.001 (-0.018) NI | +0.000 (0.000) NI |
| V5-V6 (1,913, 1,556) | +0.015 (-0.003) NI | +0.016 (-0.002) NI | +0.015 (-0.002) NI | 0.000 (0.000) NI |
| V6-V8 (2,354, 1,818) | -0.004 (-0.017) NI | -0.004 (-0.018) NI | 0.000 (-0.014) NI | 0.000 (0.000) NI |
| V7-V9c (2,626, 2,016) | +0.002 (-0.015) NI | +0.004 (-0.013) NI | -0.002 (-0.018) NI | -0.000 (-0.001) NI |

Exploratory (B10k\*). Every challenger is non-inferior on every amplicon region except bscampp at V1-V2 (bound 0.19 pp past the margin); all three fall about 2 to 3 pp short at full-length, so none meets the full-length and V4-150 headline at B10k\*. T30 extra challengers (full-length lower bound > -0.02): none qualify, so T30 runs the core arms only.

Archive: `/ddn_scratch/lxxu/phylobench-archive/t28-final/` (git `9ad3bc2`): `runs.tsv` `ee5635e57f88246601c83e7ab9710015ca832a74eb404f16d8a9ceedb1cb5d3a`, `accuracy.tsv` `d581e907008001ce125f93675f367f4cef303e48235cbe294e03852d775d7417`, `concordance.tsv` `609aeb592eebda1f063586d6faa33a1e7441a5b4573fb83b2546ea15a8ef49df`, `t28-within3.txt` `1967f34fede1f90fa0ded30c007170d458eeb7c1357b5da433e8f6f9c2c5d3f1`, `regions.tgz` `6efb4c54908afada574d7ae5df20d8a966246c8833242038784687965c13931c`. T22 tables kept as `*-t22-final.tsv` (2026-11-10).

## Two-failure rule

Same as `docs/runbook-t22.md`: first failure of a job class `(rule, backbone)` is diagnosed and resubmitted with `--rerun-incomplete`; the second stops the stage; a harness DNF (`oom`, `timeout`, `error`) is data, not a failure. If the tmux session dies, paste `squeue -u lxxu` before restarting anything.
