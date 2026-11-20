# T29 runbook: accuracy levers at B10k\* on `gg2-wol`

Stages S0 (UPP probe at B10k\*, done), S0b (UPP sizing at Bfull\*), S1 (36 cells), S2 (gated, 9 cells). Same conventions as `docs/runbook-t28.md`: run from `barnacle-login`, in `/ddn_scratch/lxxu/phylobench`, `phylobench` env active, inside `tmux`, SLURM account `knightlab`, archive directory `/ddn_scratch/lxxu/phylobench-archive`. One command per block; after each step, paste back what the "Paste" line says. Every `snakemake` run goes through a profile (`profiles/slurm-pilot`), targets come before `--configfile` (which takes several files and would swallow them), and the user runs everything; agents never touch the cluster. Never run two `snakemake` invocations in the same working directory. Cluster facts: `docs/cluster.md`.

What T29 asks: do two levers close the T28 full-length gap of `bscampp` against `sepp-prod` at B10k\* (about -3 pp within 3 edges)? Lever 1: UPP alignment instead of `hmmalign` (`upp/bscampp`). Lever 2: BSCAMPP subtree size 5000 (`hmmalign/bscampp-b5000`). Tool notes: `docs/tools/upp.md`, `docs/tools/bscampp.md`. Plan and decisions: `TICKETS.md` T29 and `design.md` section 9 item 10.

## Pre-registered rules

- **Primary table** (T30 rule): each T29 cell's within-3 difference against `raw/sepp-prod`; non-inferior iff the one-sided 95% lower bound is above -0.02 (cluster bootstrap over held-out groups, B = 10,000, seed 0; unplaced and DNF are misses; a DNF reference is undecidable).
- **Lever effect:** a lever is **shown** iff at **full-length** the paired within-3 difference (lever minus `hmmalign/bscampp`) has a one-sided 95% lower bound above 0, same bootstrap. `within_k_table.py --reference hmmalign/bscampp --margin 0` computes it; the verdict `non-inferior` is exactly "effect shown", `not-shown` otherwise. Other regions are descriptive. Minimum detectable gain about 1.5 to 2 pp (UNVERIFIED); a null is worded "not shown".
- **S2 gate:** `upp/bscampp-b5000` (9 cells) runs iff both levers are shown, after the user confirms.
- Exploratory: both levers are tested on the same queries that revealed the bscampp full-length deficit in T28.
- UPP alignment probably makes `upp/bscampp` e2e cost exceed `sepp-prod` at B10k\*: T29 measures an accuracy lever, not a cost win.

## Matrix

| Item | Value |
|---|---|
| Config | `config/t29.yaml` (S1), `+ config/t29-combo.yaml` (S2) |
| Dataset, backbone | `gg2-wol`, `B10k` only (no Bfull\* cell in T29) |
| Regions | T28's frozen 9: `full-length`, `V4-150`, `V3-V4`, `V1-V2`, `V1-V3`, `V4-V5`, `V5-V6`, `V6-V8`, `V7-V9c` |
| S1 cells (36) | reused from T28 (18): `raw/sepp-prod` x 9, `hmmalign/bscampp` x 9. New (18): `upp/bscampp` x 9, `hmmalign/bscampp-b5000` x 9 |
| S2 cells (+9) | `upp/bscampp-b5000` x 9 |
| Profile | `profiles/slurm-pilot` |

Reuse: the run_ids of the 18 reused cells are identical to T28's, their outputs are complete, and no existing rule's code, params, inputs or env changed, so Snakemake does not rerun them. `collect_tables` rewrites `results/*.tsv` with the planned cells only, so T28's other rows leave them (preflight archives T28's tables first).

Expected jobs in the S1 dry run (75 in all; it may show fewer, never more): `upp_align` 9, `upp_queries` 9, `bscampp_place` 9 (all reading `upp.aligned.fasta`), `bscampp_b5000_place` 9, `guppy_tog` 18, `evaluate` 18, `run_specs` 1, `collect_tables` 1, `all` 1. Absent: `sepp_prod`, `hmmalign`, `mask_hmmalign`, `similarity`, `extract_regions`, `gg2_wol_build`, `refprep_*`, `sepp_pinned_place`, `pplacer_scampp_place`, `epa_ng_scampp_place`, `oracle_region`, and anything for `gg2` or `Bfull`.

Cost estimate (UNVERIFIED for amplicons). UPP at `-A 100` took 381.8 s for 300 full-length queries (S0, 16 threads, EPYC 7302), so 3,000 queries is about 10 x 381.8 s = 3,818 s, about 1.1 h per cell at full-length, and at most about 9.6 h for the 9 cells if the amplicon regions cost as much (UPP aligns fewer columns there, so likely less). Place about 0.3 node-hours for `bscampp-b5000`, about 0.15 for `upp/bscampp`, `guppy_tog` and `evaluate` about 0.5. S1 total about 11 to 12 node-hours. Pre-registered stop: **30 node-hours** (`sacct`). S2 about 0.5 node-hours.

Note: the reused T28 cells ran on AMD EPYC 7302 while new cells may land on other CPUs; compare cost only within one `cpu_model` (the cost table prints it).

## S0 result (user ran on Barnacle at `b611d48`)

| Check | Result |
|---|---|
| Env | exactly one lock `MATCH` for `envs/sepp.linux-64.pin.txt`; `pplacer --version` prints `dev` (alpha23 proven by the lock) |
| U1 flags | PASS: `run_upp.py -h` lists `-s -a -t -A -m -x -o -d -p`; `-d` and `-p` need absolute paths |
| U2 hidden binaries | PASS: `main.config` paths all under the env `bin`; `~/.sepp` absent; `config.py` mentions it only in a docstring. `pasta path='run_pasta.py'` is unpinned but PASTA never runs because we always pass `-a` and `-t` |
| `-A 10` at B10k\*, 300 queries | timeout at 2,400 s (about 1,900 HMMs x 16 chunks, about 30k `hmmsearch` launches); U5 fails at `-A 10` |
| `-A 100` at B10k\*, 300 queries (FIXED) | ok, wall 381.8 s (U5 pass, <= 1,080), peak 4,617 MB cgroup (U4 pass), AMD EPYC 7302 |
| U3 | backbone width 1598, 10,000 rows; full width 1660; masked width 1598; 62 insertion columns; full minus insertion columns equals masked; backbone rows identical after `upper()`: approach A (use `upp_alignment_masked.fasta`, validated by `upp_query_rows`) |
| U7 | 300/300 queries present, no unknown names, 0 lowercase in masked queries, masked query characters `- A C G T` |
| Tiny fixture, `-A 10`, 2 threads | ok, 13.5 s, 3.2 GB, checker all True |
| U6 | UPP joins CI smoke; the first CI run checks `smoke_min_exact: 8` for `upp/bscampp` (if only that fails, UPP leaves smoke, threshold stays) |

`-A 100` is the single pinned setting everywhere.

## S0b: UPP at Bfull\* (sizing; on its own; not on T29's critical path)

Holds one `highmem` node, exclusive, for up to 3 h. `-n 1` is required: on an `--exclusive` 64-CPU node, `-c 16` without it starts 4 copies (see the result below); the row is Xeon Gold 6242 cost data for T30. The harness timeout is 8,400 s, one tenth of `upp_align`'s 84,000 s cap (the 300-query run times 10 approximates 3,000 queries). Needs `results/t29-s0/check_upp.py` and `results/t29-s0/b10k/q300.fasta` from S0.

```bash
tmux new -s t29b
```

Paste: nothing.

```bash
cd /ddn_scratch/lxxu/phylobench && conda activate phylobench && ls results/t29-s0/check_upp.py results/t29-s0/b10k/q300.fasta results/t29-s0/bfull
```

Paste: the output (expect both files and the `bfull` directory listed).

```bash
for d in .snakemake/conda/*/; do test -x "$d/bin/run_upp.py" || continue; if diff -q <(conda list -p "${d%/}" --explicit --md5 | grep '^https' | sort) <(grep '^https' envs/sepp.linux-64.pin.txt | sort) >/dev/null; then echo "MATCH ${d%/}"; else echo "other ${d%/}"; fi; done
```

Paste: the output. STOP unless exactly one `MATCH` line. Then set the variable to that path and prepare the absolute work directories (`-d` and `-p` need absolute paths):

```bash
UPP_ENV=<the MATCH path>; rm -rf results/t29-s0/bfull/out results/t29-s0/bfull/tmp; mkdir -p results/t29-s0/bfull/out results/t29-s0/bfull/tmp; echo "$PWD/results/t29-s0/bfull"
```

Paste: the printed path (expect it to start with `/ddn_scratch/lxxu/phylobench`).

```bash
srun -A knightlab -p highmem --exclusive -n 1 -c 16 --mem=0 -t 03:00:00 python -m phylobench.measure --out results/t29-s0/bfull/upp.run.json --threads 16 --timeout 8400 -- env PATH="$UPP_ENV/bin:$PATH" run_upp.py -s results/t29-s0/b10k/q300.fasta -a data/gg2-wol/Bfull/aligned-dna-sequences.fasta -t data/gg2-wol/Bfull/tree.nwk -A 100 -m dna -x 16 -o upp -d "$PWD/results/t29-s0/bfull/out" -p "$PWD/results/t29-s0/bfull/tmp" > logs/t29-s0-bfull.txt 2>&1; cat results/t29-s0/bfull/upp.run.json
```

Paste: the output.

```bash
srun -A knightlab -p short -n 1 -c 2 --mem=64G -t 03:00:00 python results/t29-s0/check_upp.py results/t29-s0/bfull/out data/gg2-wol/Bfull/aligned-dna-sequences.fasta results/t29-s0/b10k/q300.fasta
```

Paste: the output.

Decision (pre-registered):

- **B1 GO:** status `ok` (10 x wall <= 84,000 s, a conservative bound for 3,000 queries since the fixed HMM-ensemble cost is counted ten times), `peak_mem_mb` < 1,400,000 (a lower bound for 3,000 queries; T30's harness records the real peak, an OOM is a DNF row), checker values as at B10k\* (backbone rows identical, full minus insertion columns equals masked). T30 plans the upp COST cells.
- **B2 UNSIZED (time):** status `timeout`. This shows only T(300) > 2.33 h; 3,000 queries might still fit. T30 plans no upp cell by default and reports the S0b row as "unsized"; the user decides.
- **B3 UNSIZED (memory):** `oom`, or peak >= 1,400,000 MB. Same handling.
- **B4 STOP:** checker values fail at Bfull\* while they passed at B10k\*, or `cpu_model` is not Xeon Gold 6242.

### S0b result

**B2 UNSIZED (time)** (2026-11-18). One task (`-n 1`), b-cn-25, Intel Xeon Gold 6242: status `timeout` at 8,402 s, peak 31,102 MB. UPP was killed before it wrote `upp_alignment.fasta`, so the checker had no input (not a B4 STOP). T(300 queries) > 2.33 h at Bfull\*; T30 plans no upp cell and reports this row as the upp Bfull\* evidence.

Incident, void run: the first attempt used `--exclusive -c 16` without `-n 1`. On the exclusive 64-CPU node srun started 4 tasks (`sacct` step `2316728.0`: NTasks 4, AllocCPUS 64), so 4 UPP copies ran side by side (the log has 4 `All options` lines and 4 temp directories; timeout at 8,402 s, 30 GB). Discarded. The S0 B10k\* and tiny probes ran on shared `short` nodes as one task (step `2316714.0`: NTasks 1) and stand. The block above now has `-n 1`.

## Preflight (S1)

```bash
tmux new -s t29
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
git checkout --detach <T29_MERGE_SHA>
```

Paste: the last line of output.

```bash
git rev-parse HEAD
```

Paste: the SHA.

```bash
git diff 9ad3bc2 HEAD --stat -- envs/
```

Paste: the output (expect empty; otherwise finished jobs rerun).

```bash
python -c "import phylobench, snakemake; print(phylobench.__file__, snakemake.__version__)"
```

Paste: the output (expect a path under `/ddn_scratch/lxxu/phylobench/src` and `9.27.0`).

```bash
python -m pytest -q tests/test_upp_queries.py tests/test_t29_config.py tests/test_within_k.py tests/test_within_k_table.py tests/test_tables.py tests/test_profile.py
```

Paste: the last line.

The 18 reused cells must be complete.

```bash
for r in full-length V4-150 V3-V4 V1-V2 V1-V3 V4-V5 V5-V6 V6-V8 V7-V9c; do ls results/gg2-wol/B10k/$r/raw/sepp-prod/accuracy.tsv results/gg2-wol/B10k/$r/hmmalign/bscampp/accuracy.tsv; done 2>&1 | grep -c "accuracy.tsv$"
```

Paste: the number (expect 18; a lower count means a missing file, paste the `ls` errors).

T28's archived tables must exist before `collect_tables` overwrites `results/*.tsv`. The block prints `t28-final present` or `STOP`.

```bash
mkdir -p logs; test -s /ddn_scratch/lxxu/phylobench-archive/t28-final/accuracy.tsv && test -s /ddn_scratch/lxxu/phylobench-archive/t28-final/runs.tsv && echo "t28-final present" || echo "STOP: t28-final missing from the archive; do not continue"
```

Paste: the output. On `STOP`, do not continue.

`cp -n` never overwrites, so re-running cannot replace the T28 copy with a T29 table.

```bash
for t in runs accuracy concordance; do cp -n results/$t.tsv /ddn_scratch/lxxu/phylobench-archive/$t-t28-final.tsv; done; ls -l /ddn_scratch/lxxu/phylobench-archive/*-t28-final.tsv; sha256sum /ddn_scratch/lxxu/phylobench-archive/*-t28-final.tsv
```

Paste: the output.

## S1: B10k\* (36 cells; 18 new)

Dry run to a file.

```bash
snakemake all -n --profile profiles/slurm-pilot --configfile config/t29.yaml > logs/t29-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A30 "^Job stats" logs/t29-dryrun.txt | head -40
```

Paste: the output. STOP if it shows any job for `sepp_prod`, `hmmalign`, `mask_hmmalign`, `similarity`, `extract_regions`, `gg2_wol_build`, `refprep_*`, or another place rule, or more than the expected counts (Matrix section).

```bash
snakemake all --list-changes code --profile profiles/slurm-pilot --configfile config/t29.yaml; echo "--- params"; snakemake all --list-changes params --profile profiles/slurm-pilot --configfile config/t29.yaml; echo "--- input"; snakemake all --list-changes input --profile profiles/slurm-pilot --configfile config/t29.yaml
```

Paste: the output. Expect `code` empty, `params` only `run_specs`, `input` only `collect_tables`; anything else is a STOP (a reused cell would rerun).

Wait for my confirmation, then record the start time.

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t29.start
```

Paste: nothing.

```bash
snakemake all --profile profiles/slurm-pilot --jobs 6 --keep-going --configfile config/t29.yaml 2>&1 | tee logs/t29.log
```

Paste: nothing until it ends; then the last 30 lines of `logs/t29.log`.

```bash
grep -c "Error in rule" logs/t29.log
```

Paste: the number, and the matching lines if it is not 0. Two-failure rule: if two jobs of one class `(rule, backbone)` fail, stop the stage and paste; do not rerun. A UPP job past 84,000 s is a DNF row (data), not a failure.

Stage gate.

```bash
python workflow/scripts/check_stage.py config/t29.yaml
```

Paste: the full output. `FAIL` lines block the stage; `INFO` lines are OOM or timeout DNFs (data).

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t29.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number (stop and paste above 30).

## Reuse check and within-3 table

The 18 reused cells must be identical to T28's archived rows.

```bash
python - <<'EOF'
import pandas as pd
new = pd.read_csv("results/accuracy.tsv", sep="\t")
old = pd.read_csv("/ddn_scratch/lxxu/phylobench-archive/t28-final/accuracy.tsv", sep="\t")
ids = sorted(r for r in set(new["run_id"]) if r.startswith("gg2-wol/B10k/") and r.endswith(("/raw/sepp-prod", "/hmmalign/bscampp")))
print("reused run_ids", len(ids))
key = ["run_id", "query"]
a = new[new["run_id"].isin(ids)].sort_values(key).reset_index(drop=True)
b = old[old["run_id"].isin(ids)].sort_values(key).reset_index(drop=True)
print("rows", len(a), len(b), "identical:", a.equals(b))
EOF
```

Paste: the output (expect 18 run_ids and `identical: True`; otherwise STOP).

Primary table: every T29 cell against `raw/sepp-prod`.

```bash
python workflow/scripts/within_k_table.py config/t29.yaml | tee results/t29-within3.txt
```

Paste: the full output (27 rows: 36 cells less the 9 `sepp-prod` references). The `hmmalign/bscampp` rows must equal T28's, e.g. full-length `-0.030 (-0.052)`.

## Lever table and gate

Each lever against `hmmalign/bscampp`, margin 0.

```bash
python workflow/scripts/within_k_table.py config/t29.yaml --reference hmmalign/bscampp --margin 0 | tee results/t29-levers.txt
```

Paste: the full output (18 rows: `upp/bscampp` and `hmmalign/bscampp-b5000` x 9 regions).

```bash
awk 'NR==1 || $3=="full-length"' results/t29-levers.txt; awk 'NR>1 && $3=="full-length" && $NF=="non-inferior"{n++} END{print n+0, "of 2 levers shown at full-length; S2 runs iff 2"}' results/t29-levers.txt
```

Paste: the output. Wait for my confirmation before S2. Wording: `non-inferior` at margin 0 means "effect shown" (lower bound > 0); `not-shown` means "not shown", never "no effect".

## Cost table

End-to-end cost of the new cells beside `hmmalign/bscampp`, with the one-time `hmmbuild` row (hmmalign's per-backbone model cost, excluded from e2e, against UPP's per-region ensemble inside its align row). Compare only within one `cpu_model`.

```bash
python - <<'EOF' | tee results/t29-cost.txt
import pandas as pd
from phylobench.summaries import end_to_end
runs = pd.read_csv("results/runs.tsv", sep="\t")
runs = runs[(runs["dataset"] == "gg2-wol") & (runs["backbone"] == "B10k")]
cols = ["region", "aln", "arm", "query_n", "place_s", "e2e_s", "peak_mem_mb", "cpu_model", "e2e_status"]
print(end_to_end(runs)[cols].sort_values(["region", "arm", "aln"]).to_string(index=False))
align = runs[runs["arm"].isna() & runs["region"].notna() & (runs["phase"] == "align")]
print(align[["region", "aln", "wall_s", "peak_mem_mb", "cpu_model", "status"]].sort_values(["aln", "region"]).to_string(index=False))
print(runs[runs["arm"].isna() & runs["region"].isna() & (runs["phase"] == "align")][["tool", "wall_s", "peak_mem_mb", "cpu_model", "status"]].to_string(index=False))
EOF
```

Paste: the full output (first table: e2e per cell; second: align rows per region, `hmmalign` and `upp`; third: the `hmmbuild` row).

## S2 (only if the gate passes and I confirm)

```bash
snakemake all -n --profile profiles/slurm-pilot --configfile config/t29.yaml config/t29-combo.yaml > logs/t29-s2-dryrun.txt 2>&1; grep -A20 "^Job stats" logs/t29-s2-dryrun.txt
```

Paste: the output. Expect `bscampp_b5000_place` 9, `guppy_tog` 9, `evaluate` 9, `run_specs`/`collect_tables`/`all` 1 each; no `upp_align` or `upp_queries`. Anything else is a STOP.

```bash
snakemake all --profile profiles/slurm-pilot --jobs 6 --keep-going --configfile config/t29.yaml config/t29-combo.yaml 2>&1 | tee logs/t29-s2.log; grep -c "Error in rule" logs/t29-s2.log
```

Paste: the last 30 lines and the error count.

```bash
python workflow/scripts/check_stage.py config/t29.yaml config/t29-combo.yaml
```

Paste: the full output.

```bash
python workflow/scripts/within_k_table.py config/t29.yaml config/t29-combo.yaml | tee results/t29-within3-s2.txt; python workflow/scripts/within_k_table.py config/t29.yaml config/t29-combo.yaml --reference hmmalign/bscampp --margin 0 | tee results/t29-levers-s2.txt
```

Paste: the full output of both.

## Archive and copy to the laptop

```bash
git rev-parse HEAD > results/GIT_SHA
```

Paste: nothing.

```bash
mkdir -p /ddn_scratch/lxxu/phylobench-archive/t29-final
```

Paste: nothing.

```bash
cp results/runs.tsv results/accuracy.tsv results/concordance.tsv results/t29-*.txt results/GIT_SHA /ddn_scratch/lxxu/phylobench-archive/t29-final/
```

Paste: nothing.

```bash
tar czf /ddn_scratch/lxxu/phylobench-archive/t29-final/logs.tgz config/t29.yaml config/t29-combo.yaml logs/t29* logs/t29-s0-*; echo $?
```

Paste: the number.

```bash
cd /ddn_scratch/lxxu/phylobench-archive/t29-final && sha256sum runs.tsv accuracy.tsv concordance.tsv t29-*.txt logs.tgz > SHA256SUMS && cat SHA256SUMS
```

Paste: the output.

On the laptop, from the repo root. The first command must print nothing.

```bash
test -e results/t29-final && echo "STOP: results/t29-final already exists"
```

Paste: the output.

```bash
scp -r lxxu@barnacle.ucsd.edu:/ddn_scratch/lxxu/phylobench-archive/t29-final results/ && cd results/t29-final && shasum -a 256 -c SHA256SUMS
```

Paste: the output.

## Result (2026-11-19, at `ebc5966`)

S1: 75 of 75 jobs, 0 `Error in rule`, `check_stage.py` all PASS (no DNFs), 2.8 node-hours (stop 30), 13.7 h wall (mostly queue). Reuse check: 18 reused run_ids, 43,874 rows identical to `t28-final`. Archive `t29-final/` (tables, `t29-*.txt`, `logs.tgz`, `SHA256SUMS`, `GIT_SHA`), checksums verified on the laptop. Every cell ran on AMD EPYC 7302.

Exploratory: both levers were tested on the same queries that revealed the T28 bscampp full-length deficit.

**Lever table** (`results/t29-levers.txt`; lever minus `hmmalign/bscampp`, fraction within 3 edges; lower = one-sided 95% bound, cluster bootstrap B = 10,000; shown iff lower > 0):

| Region (n) | `hmmalign/bscampp-b5000` | `upp/bscampp` |
|---|---|---|
| **full-length (3,000)** | **+0.016 (+0.0068) shown** | +0.0047 (-0.0007) not shown |
| V4-150 (2,703) | +0.018 (+0.0052) | +0.0015 (-0.0061) |
| V3-V4 (2,542) | -0.0098 (-0.022) | +0.0020 (-0.0034) |
| V1-V2 (2,119) | +0.015 (+0.0070) | +0.014 (+0.0032) |
| V1-V3 (2,058) | -0.0063 (-0.016) | +0.0092 (0.0000) |
| V4-V5 (2,622) | -0.0050 (-0.019) | +0.0015 (-0.0030) |
| V5-V6 (1,913) | +0.0021 (-0.0067) | -0.0068 (-0.025) |
| V6-V8 (2,354) | +0.0051 (-0.0013) | +0.0064 (-0.0017) |
| V7-V9c (2,626) | -0.015 (-0.026) | +0.0065 (+0.0007) |

Only full-length decides; other regions are descriptive. **S2 gate: 1 of 2 levers shown, so S2 did not run.** "Not shown" for UPP means not shown, not "no effect"; its full-length gain is under about 1 pp.

**Primary table** (`results/t29-within3.txt`; minus `raw/sepp-prod`, non-inferior iff lower > -0.02). Full-length: `hmmalign/bscampp` -0.030 (-0.052), equal to T28; `hmmalign/bscampp-b5000` -0.014 (-0.036); `upp/bscampp` -0.025 (-0.046); none non-inferior. Every amplicon region is non-inferior for all three except `hmmalign/bscampp` V1-V2 (-0.022) and `hmmalign/bscampp-b5000` V3-V4 (-0.028), V4-V5 (-0.022), V7-V9c (-0.030).

**Cost** (`results/t29-cost.txt`, EPYC 7302, 16 threads; e2e = align + place + convert, `hmmbuild` 3.4 s excluded):

| Full-length cell | e2e (s) | place (s) | peak (MB) |
|---|---|---|---|
| `hmmalign/bscampp` | 115 | 52 | 20,664 |
| `hmmalign/bscampp-b5000` | 138 | 75 | 54,476 |
| `upp/bscampp` | 3,248 | 62 | 20,399 |
| `raw/sepp-prod` | 741 | (one phase) | 98,257 |

UPP align: 3,186 s and 4,694 MB at full-length, 312 to 835 s and 1.8 to 2.4 GB per amplicon region, against `hmmalign` 5 to 62 s. `upp/bscampp` costs 4.4x `sepp-prod` end to end at full-length for a gain that is not shown.

For context, `sepp-prod`'s own log at full-length (T28 run, 735 s cumulative): HMM ensemble build, search and align 20:27:24 to 20:30:27 (about 183 s, 25%), pplacer to 20:39:12 (about 525 s, 71%), merge 27 s. SEPP's time is mostly placement, consistent with the subtree-size lever mattering more than the alignment lever.

What it means for T30: nothing changes. Both levers stay COST cells; no lever is an extra challenger; with S0b unsized and S2 not run, T30 has no upp cell.
