# T31 runbook: WoL genome held-out set (`gg2-wol`) on Barnacle

Stages S0 (read-only probe, before the T31 PR merges), build, ref-prep, checksum. Same conventions as `docs/runbook-t22.md`: run from `barnacle-login`, in `/ddn_scratch/lxxu/phylobench`, `phylobench` env active, inside `tmux`, SLURM account `knightlab`, archive directory `/ddn_scratch/lxxu/phylobench-archive`. One command per block; after each step, paste back what the "Paste" line says. Never run two `snakemake` invocations in the same working directory. Cluster facts: `docs/cluster.md`.

What the set is and how truth is worded: `docs/gg2.md` ("WoL genome held-out set (T31)"). Decision rule it feeds: `design.md` section 6 and `TICKETS.md` T30.

Prerequisite: `data/gg2/` exists (T22): `data/gg2/shipped/` and the `Bfull`/`B10k` backbones.

## Pre-registered stop/go (fixed before any gg2-wol data)

| Rule | Read from | Decision |
|---|---|---|
| G1 | P1 (G leaves) | G >= 10,000: N = 3,000. 6,000 to 9,999: N = 2,000. Below 6,000: STOP |
| G2 | P2 | One leaf per genome, else STOP (default policy: hold out all copies, query one) |
| G3 | P9 | STOP if clade queries < 500, or G singletons that can be drawn outside the clades < N minus clade queries |
| G4 | P8 | STOP if V4-150 primer hits on G would give < 2,000 queries at N (power) |
| G5 | P8 | A `c` region with fewer than 30 primer-hit G rows is dropped |
| G6 | `sinfo -p highmem` | highmem node CPU models decide routing: `--partition` if all are Xeon 6242, else `--nodelist` |

`probe_wol.py` prints G1 to G5 verdicts from the numbers; G2 prints `REVIEW` when any G-like label fails the strict `^G\d{9}$` pattern (copy suffixes), and a person reads P2 and P3. P9 uses the default clade counts in `gg2.WOL_RANKS` (phylum 6, class 9, order 12, family 18, genus 18); if it reports a shortfall at a rank, S0 sets the final counts (a one-line change in the PR) before the build.

## S0: probe (read-only, no benchmark jobs)

```bash
tmux new -s t31
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
git checkout --detach <T31_PR_SHA>
```

Paste: the last line.

```bash
git rev-parse HEAD
```

Paste: the SHA.

```bash
python -m pytest -q tests/test_wol.py tests/test_heldout.py tests/test_gg2.py
```

Paste: the last line.

```bash
mkdir -p logs
```

Paste: nothing.

The probe holds the 331k-row alignment (about 0.5 GB) and a 16 * R * L byte similarity table (about 8.5 GB); 96 GB is generous.

```bash
srun -A knightlab -p short -c 8 --mem=96G -t 04:00:00 python workflow/scripts/probe_wol.py --reference-dir data/gg2/shipped/sepp-reference --taxonomy data/gg2/shipped/tax/taxonomy.tsv --backbone-dirs data/gg2/Bfull data/gg2/B10k 2>&1 | tee logs/t31-s0-probe.txt
```

Paste: the full output (P1 to P10 and the G1 to G5 verdict lines).

```bash
sinfo -p highmem -N -o "%N %T %m %c"
```

Paste: the output.

`srun -p highmem` lands on one arbitrary node, so ask each node by name (a node that is busy queues until it is free).

```bash
for n in $(sinfo -p highmem -N -h -o "%N" | sort -u); do echo "$n: $(srun -A knightlab -p highmem -w $n -n 1 -t 5 lscpu | grep 'Model name')"; done
```

Paste: the output (one line per node). G6: if every highmem node is a Xeon 6242, route by partition; otherwise tell me and the T30 jobs pin a node.

Stop here and wait for my go/stop on G1 to G6. No build before that, and the PR does not merge before a go.

### S0 result (2026-11-06, `41b7ce6`): GO

| Rule | Probe | Decision |
|---|---|---|
| G1 | 12,283 G leaves of 331,269 (current Bfull\* 11,743, B10k\* 344) | GO, N = 3,000 |
| G2 | `REVIEW` listed only `RS_`/`GB_` record names whose contig accession contains `G` and nine digits (e.g. `..._NZ_JAGYPG#.1`); none is a `G` leaf with a copy suffix; P3: 0 G leaves with an `RS_`/`GB_` twin | GO, one leaf per genome |
| G3 | P9 with the default `WOL_RANKS`: 63 clades, 805 clade queries, 2,195 singletons wanted of 10,506 G outside clades, draw ok | GO, `WOL_RANKS` unchanged (P6 has at least 49/103/221/318/278 clades per rank) |
| G4 | V4-150 primer hits on 10,620 of 12,283 G rows | GO, about 2,594 V4-150 queries expected |
| G5 | Every `c` window has thousands of rows | No `c` region dropped |
| G6 | `b-cn-25`, `b-cn-26`: Intel Xeon Gold 6242 | Route by `--partition highmem` |

Expected queries per primer region at N = 3,000 (primer-hit rate on all G rows): V4-150 2,594, V3-V4 2,453, V4-V5 2,551, V6-V8 2,319, V1-V2 2,019, V1-V3 1,962, V5-V6 1,906, V7-V9 99 (1492R sits at the 3' end, which genome-derived genes often lack; `V7-V9c` keeps about 2,614). Descriptive: G length median 1,497 (min 769, q10 1,140, q90 1,531); 226 of 2,000 sampled G leaves have a >=99.9 neighbour (seeding target 100); lineage depth 7 for 11,606 of 12,283.

## Build (after the PR merged and S0 said go)

```bash
git fetch origin && git checkout --detach <T31_MERGE_SHA> && git rev-parse HEAD
```

Paste: the SHA.

Dry run to a file (targets go before `--configfile`, which would otherwise swallow them); the build is the only new job (`gg2_download` and `gg2_unzip` are already done).

```bash
snakemake data/gg2-wol/MANIFEST.tsv -n --profile profiles/slurm --configfile config/full.yaml > logs/t31-build-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A12 "^Job stats" logs/t31-build-dryrun.txt
```

Paste: the output. Expect exactly one job, `gg2_wol_build`; none for `gg2_build`, `gg2_download`, `gg2_unzip`.

```bash
snakemake data/gg2-wol/MANIFEST.tsv --profile profiles/slurm --configfile config/full.yaml 2>&1 | tee logs/t31-build.log
```

Paste: the last 15 lines. Expected wall: under an hour (T18's build took about the same).

```bash
tail -n 12 logs/gg2-wol/build.log
```

Paste: the output (held-out leaf count, queries, seeded, duplicates removed, per-region bin lines, any `WARNING`).

```bash
python -c "
import pandas as pd
h = pd.read_csv('data/gg2-wol/heldout.tsv', sep='\t')
q = h[h['query'] == 1]
print(len(q), 'queries;', q['group_type'].value_counts().to_dict())
print('non-G queries:', (~q['leaf'].str.fullmatch(r'G\d{9}')).sum())
print('clades:', h[h.group_type == 'clade'].groupby('rank')['group'].nunique().to_dict())
print('seeded:', int(h['seeded'].sum()))
"
```

Paste: the output. `non-G queries` must be 0.

```bash
awk -F'\t' '$1=="bin" && ($3=="full-length" || $3=="V4-150") {print $2, $3, $4, $6, $8}' data/gg2-wol/MANIFEST.tsv
```

Paste: the output (bin counts per backbone for full-length and V4-150; flagged bins are data, not errors).

## Ref-prep (FastTree me-tree, epa-ng model, hmmbuild; both backbones)

Per-backbone prep runs fresh for gg2-wol (the gg2 prep is for a different backbone). Bfull FastTree used 196 GB in T22, so the profile's `long` resources apply.

```bash
snakemake results/ref-prep/gg2-wol/B10k/me-tree.nwk results/ref-prep/gg2-wol/B10k/epa-ng.model results/gg2-wol/B10k/align/backbone.hmm results/ref-prep/gg2-wol/Bfull/me-tree.nwk results/ref-prep/gg2-wol/Bfull/epa-ng.model results/gg2-wol/Bfull/align/backbone.hmm -n --profile profiles/slurm --configfile config/full.yaml > logs/t31-prep-dryrun.txt 2>&1
```

Paste: nothing.

```bash
grep -A14 "^Job stats" logs/t31-prep-dryrun.txt
```

Paste: the output. Expect per backbone one each of `refprep_me_tree`, `refprep_epa_ng_model`, `backbone_stockholm`, `hmmbuild`, and no `gg2_*`.

Wait for my confirmation, then record the start time.

```bash
date +%Y-%m-%dT%H:%M:%S > logs/t31-prep.start
```

Paste: nothing.

```bash
snakemake results/ref-prep/gg2-wol/B10k/me-tree.nwk results/ref-prep/gg2-wol/B10k/epa-ng.model results/gg2-wol/B10k/align/backbone.hmm results/ref-prep/gg2-wol/Bfull/me-tree.nwk results/ref-prep/gg2-wol/Bfull/epa-ng.model results/gg2-wol/Bfull/align/backbone.hmm --profile profiles/slurm --jobs 8 --keep-going --configfile config/full.yaml 2>&1 | tee logs/t31-prep.log
```

Paste: the last 30 lines.

```bash
grep -c "Error in rule" logs/t31-prep.log; grep -h '"status"' results/ref-prep/gg2-wol/*/me-tree.run.json results/gg2-wol/*/align/hmmbuild.run.json
```

Paste: the output (expect 0 errors and four `"status": "ok"`).

```bash
sacct -X -n -P -u lxxu -A knightlab -S "$(cat logs/t31-prep.start)" --format=Elapsed,NNodes | awk -F'|' '{n=split($1,d,"-"); days=(n==2)?d[1]:0; t=(n==2)?d[2]:d[1]; split(t,h,":"); s+=(days*24+h[1]+h[2]/60+h[3]/3600)*$2} END{printf "%.1f node-hours\n", s}'
```

Paste: the number (estimate 5 to 10 node-hours for build plus prep, UNVERIFIED).

## Checksum and archive

```bash
sha256sum data/gg2-wol/MANIFEST.tsv data/gg2-wol/heldout.tsv data/gg2-wol/queries.fasta data/gg2-wol/Bfull/tree.nwk data/gg2-wol/B10k/tree.nwk
```

Paste: the output (I record these in `docs/gg2.md`).

```bash
tar czf /ddn_scratch/lxxu/phylobench-archive/t31-gg2-wol.tgz data/gg2-wol logs/t31-*.txt logs/t31-*.log logs/gg2-wol && ls -l /ddn_scratch/lxxu/phylobench-archive/t31-gg2-wol.tgz
```

Paste: the output.

## What to paste back, in order

1. S0: SHA, pytest line, full probe output, `sinfo` and `lscpu` lines. I decide G1 to G6.
2. Build: dry-run job stats, last log lines, the held-out summary, the bin counts.
3. Ref-prep: dry-run job stats, error count and statuses, node-hours.
4. Checksums and the archive listing.
