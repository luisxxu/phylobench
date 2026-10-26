# T21 HPC profile and pilot

Pilot of every arm at B10k\* on Barnacle (SLURM, `knightlab`), plus a Bfull\* DIPPER GPU probe. Its output is the Bfull\* go/no-go table that fixes `config/full.yaml` (T22). Cluster facts are in `docs/cluster.md`; the measurement contract is `docs/harness.md`.

## Profile (`profiles/slurm/config.yaml`)

- Executor `slurm`, `retries: 0` (a DNF is data), `use-conda: true`, `jobs: 20`, no `group` keys. One rule per sbatch gives each job its own `srun` step cgroup, so `mem_source=cgroup` is primary (`docs/harness.md` section 4).
- Defaults (`evaluate`, `collect_tables`, data prep): `short`, 1440 min, 32 GB, non-exclusive. `evaluate` therefore runs beside other jobs and is excluded from all timing.
- The 10 CPU place rules (`sepp_prod` and every `*_place`): `long`, 4320 min (72 h), 500000 MB, `--exclusive`. A 16-thread placement gets a whole node, so place timings are free of neighbours.
- Prep and convert rules (`hmmbuild`, `hmmalign`, `refprep_me_tree`, `guppy_tog`): `long`, 4320 min, 196608 MB, shared nodes. Barnacle pins each job to its cores (`TaskPlugin=task/cgroup,task/affinity`, `CR_CORE_MEMORY`), but neighbours still share memory bandwidth, cache and Lustre, so these phases' wall times carry a shared-node caveat; memory stays exact (own cgroup). Changed mid-pilot (2026-10-22): with every prep and convert job waiting for a whole node, the r1 queue reached days. Stage 1's B10k\* `hmmbuild`/`hmmalign` rows ran exclusive. `sepp_prod` is one rule, so its align and convert stay exclusive while other arms' do not; end-to-end comparisons against sepp-prod carry this caveat.
- `dipper_m_gpu_place`: `gpu` partition, `gpu: 1`, `--exclusive --nodelist=b-gpu-21`, 4320 min, 196608 MB. The `gpu` partition sets `MaxMemPerCPU=12288`, so 16 CPUs cap the request at 16 x 12288 MB; asking for more makes SLURM raise the CPU count past the node's 32 and reject the job (first stage-2 submission, 2026-10-22). The plugin's `gpu` resource is used instead of `gres` because `gres` is not propagated to the jobstep `srun`.
- `profiles/slurm-pilot/config.yaml` (stage 4 and repeats r2/r3): identical except every `long` rule (exclusive place rules and the shared prep/convert rules) runs on `short` with 1440 min; on 2026-10-23 even the shared `long` jobs were estimated four days out. `short` has PriorityTier 3 and `long` 2 over the same nodes, so on a busy cluster an exclusive `long` job waited about 2 days while the same job on `short` was estimated within a day (2026-10-22, `sbatch --test-only`). The harness timeout stays 71 h (changing `timeout_s` changes every rule's shell and would rerun stages 1-2), so a B10k\* job over 24 h would be killed by SLURM without a `run.json`; such a cell is recorded here as DNF by hand. The longest B10k\* place cell in stages 1-2 took 6 min. `--jobs 6` limits how many whole nodes the pilot holds at once.
- `dipper-m.in.fasta` is `temp()` (`workflow/rules/dipper.smk`), so once one DIPPER arm has run alone it is deleted; a later run that needs it for the other `-i m` arm rebuilds it and Snakemake reruns the finished arm too (stage 4 reran `dipper-m-gpu`, ~1 GPU-minute; stage 2 rows are kept in `pilot-s2.tgz`). T22 should run `dipper-m` and `dipper-m-gpu` in the same invocation.
- `--no-requeue` is passed through `slurm_extra` rather than the plugin's `slurm-no-requeue` setting, so the profile has no plugin-specific keys and `snakemake --executor local --profile profiles/slurm` works without the plugin (used by `tests/test_profile.py`).
- `envs/phylobench.yaml` pins `snakemake-executor-plugin-slurm=2.8.0` and `-jobstep=0.6.1`. Job-side Snakemake runs under the submitting env's Python (the plugin uses `sys.executable` and `--export=ALL`), so submit from the activated `phylobench` env.

## Choices

- **Card:** RTX 3090 on `b-gpu-21` only (24 GB, 32 CPUs, 503 GiB; design section 5). GPU runtimes are comparable only within one card, so there is no fallback to another node.
- **`--force-binary`** on `dipper-m-gpu` only: the image otherwise collapses near-zero internal edges into multifurcations, which changes topology relative to the backbone and breaks the adapter's RF check. The CPU 0.1.5 build always collapses exact-zero internal edges and has no flag.
- **`timeout_s`** (`config/schema.yaml`, default 255600 s = 71 h): appended to every harness call so the harness records a DNF before the 72 h SLURM limit. A rule's own later `--timeout` wins (pplacer 86400; smoke 1).

## Configs

| Config | Cells |
|---|---|
| `config/pilot.yaml` | gg2 and emp-10k, B10k, full-length and V4-150 (emp: V4-150 only), hmmalign: 33 cells, all 11 arms |
| `config/pilot-cv.yaml` | gg2 / B10k / V4-150 / hmmalign: 11 cells, for timing repeats |
| `config/pilot-bfull-probe.yaml` | gg2 / Bfull / V4-150 / hmmalign, run with `--config arm=dipper-m-gpu` |

## Staging

Run from the activated `phylobench` env, in order. Do not start a stage until the previous stage's gates pass. The profile already sets `use-conda`, so `--use-conda` is not repeated. The archive directory is created first, and `results/runs.tsv` is regenerated on every invocation, so archive before each forced run.

```sh
mkdir -p /ddn_scratch/lxxu/phylobench-archive
```

Dry run:

```sh
snakemake all -n --profile profiles/slurm --configfile config/pilot.yaml
```

Stage 1:

```sh
snakemake all --profile profiles/slurm --config arm=sepp-prod --configfile config/pilot.yaml
```

```sh
tar czf /ddn_scratch/lxxu/phylobench-archive/pilot-s1.tgz results
```

Stage 2:

```sh
snakemake all --profile profiles/slurm --config arm=dipper-m-gpu --configfile config/pilot.yaml
```

```sh
tar czf /ddn_scratch/lxxu/phylobench-archive/pilot-s2.tgz results
```

Stage 3:

```sh
snakemake all --profile profiles/slurm --config arm=dipper-m --configfile config/pilot.yaml
```

```sh
tar czf /ddn_scratch/lxxu/phylobench-archive/pilot-s3.tgz results
```

Stage 4 (all arms, timing run r1):

```sh
snakemake all --profile profiles/slurm-pilot --jobs 6 --configfile config/pilot.yaml
```

```sh
tar czf /ddn_scratch/lxxu/phylobench-archive/pilot-r1.tgz results
```

```sh
cp results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-r1.tsv
```

Timing repeat r2 (the 11 place rules on the CV config):

```sh
snakemake all --profile profiles/slurm-pilot --jobs 6 --forcerun sepp_prod pplacer_place epa_ng_place bscampp_place epa_ng_scampp_place pplacer_scampp_place apples2_place dipper_m_place dipper_r_place dipper_m_gpu_place sepp_pinned_place --configfile config/pilot-cv.yaml
```

```sh
tar czf /ddn_scratch/lxxu/phylobench-archive/pilot-r2.tgz results
```

```sh
cp results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-r2.tsv
```

Timing repeat r3 (the 11 place rules on the CV config):

```sh
snakemake all --profile profiles/slurm-pilot --jobs 6 --forcerun sepp_prod pplacer_place epa_ng_place bscampp_place epa_ng_scampp_place pplacer_scampp_place apples2_place dipper_m_place dipper_r_place dipper_m_gpu_place sepp_pinned_place --configfile config/pilot-cv.yaml
```

```sh
tar czf /ddn_scratch/lxxu/phylobench-archive/pilot-r3.tgz results
```

```sh
cp results/runs.tsv /ddn_scratch/lxxu/phylobench-archive/runs-r3.tsv
```

Bfull\* GPU probe (Option A; Option B is off):

```sh
snakemake all --profile profiles/slurm --config arm=dipper-m-gpu --configfile config/pilot-bfull-probe.yaml
```

```sh
tar czf /ddn_scratch/lxxu/phylobench-archive/pilot-bfull-gpu-probe.tgz results
```

Summary:

```sh
python workflow/scripts/pilot_summary.py /ddn_scratch/lxxu/phylobench-archive/runs-r1.tsv /ddn_scratch/lxxu/phylobench-archive/runs-r2.tsv /ddn_scratch/lxxu/phylobench-archive/runs-r3.tsv --out results/pilot-timing.tsv
```

Columns: `dataset backbone region aln arm n n_ok wall_mean wall_cv peak_mem_max`, place phase only; CV is the sample SD over ok rows divided by the mean, NaN with fewer than two.

## Stage gates

Every stage must show, in `results/runs.tsv` and SLURM accounting:

- all rows `status=ok`
- `mem_source=cgroup`
- the step line in `sacct` shows `AllocCPUS` 16
- `cpu_s / wall_s <= 16.8` (16 threads plus 5% slack)
- GPU arm additionally: `node` is `b-gpu-21`, `gpu_peak_mb > 0`, and `dipper_to_jplace` passes its RF check

Zero-length internal edge pre-flight: count zero-length internal edges in the B10k\* backbone trees before stage 3. If nonzero and CPU DIPPER fails the RF check, report it as a finding for a user decision. Do not patch the trees or the adapter silently.

## Go/no-go criteria (pre-registered)

- pplacer at Bfull\* runs only as the capped attempt in T22 (`config/pplacer-bfull-cap.yaml`). `epa-ng` and `dipper-r` run at B10k\* only.
- Projection to Bfull\* assumes cost linear in backbone size (about 32x from B10k\*), with a per-arm query-count slope fitted from gg2 V4-150 (1k queries) and emp-10k (10k queries) at B10k\*. This is an assumption, labelled as such in the table.
- **GO** for an arm iff all its pilot rows are ok, `cpu_s / wall_s <= 16.8`, projected Bfull\* place wall for emp-100k is at most 36 h, and projected peak memory is at most 0.8 x node memory. The GPU arm additionally needs the Bfull\* probe's measured `gpu_peak_mb` at most 0.8 x 24 GB, else "unverified".
- Projected wall between 36 h and 71 h: "go, capped". Above 71 h: no-go.
- `apples2` results always carry the uDance caveat (GG2 truth was seeded by APPLES-2).
- Dropping any arm from `config/full.yaml` needs user approval in T22.

## Budget

About 60 node-hours by default; ceiling 150. Record actual node-hours below.

## Results

Barnacle, 2026-10-22 to 2026-10-26. All 33 B10k\* place cells `ok` (acceptance met); every place row `mem_source=cgroup`; max `cpu_s / wall_s` 15.7. Archives: `/ddn_scratch/lxxu/phylobench-archive/` (`pilot-s1`, `pilot-s2`, `pilot-r1..r3` `.tgz`, `runs-r1..r3.tsv`).

### Per-arm runtime and memory at B10k\*

Place phase (SEPP arms: their combined `align+place[+convert]` phase), r1; gg2 V4-150 is the mean of 3 runs.

| arm | gg2 full-length wall (s) / peak (MB) | gg2 V4-150 | emp-10k V4-150 |
|---|---|---|---|
| sepp-prod | 337.6 / 102,337 | 68.9 / 11,613 | 242.6 / 12,212 |
| sepp-pinned | 286.9 / 103,229 | 73.1 / 11,551 | 207.1 / 12,133 |
| pplacer | 248.3 / 14,795 | 26.4 / 2,500 | 239.1 / 3,546 |
| pplacer-scampp | 934 / 41,414 | 188.1 / 7,508 | 1153.3 / 25,535 |
| epa-ng-scampp | 288.8 / 37,729 | 102.4 / 19,111 | 613.6 / 33,811 |
| bscampp | 30.7 / 20,943 | 8.35 / 9,267 | 37.4 / 9,671 |
| epa-ng | 23.1 / 15,811 | 10.5 / 6,220 | 22.8 / 6,290 |
| apples2 | 26.1 / 737 | 7.27 / 603 | 42.1 / 885 |
| dipper-m | 4.36 / 194 | 3.14 / 194 | 44.3 / 252 |
| dipper-m-gpu | 4.78 / 585 | 4.27 / 590 | 11 / 631 |
| dipper-r | 6.51 / 340 | 4.45 / 340 | 69.5 / 481 |

`dipper-m-gpu` ran on `b-gpu-21` (RTX 3090) with `gpu_peak_mb` 336 (gg2) and 348 (emp-10k).

### Timing CV (3 runs, gg2 V4-150)

| arm | n_ok | wall_mean (s) | wall_cv |
|---|---|---|---|
| sepp-prod | 3 | 68.9 | 11.3% (flag: >10%) |
| sepp-pinned | 3 | 73.1 | 17.1% (flag: >10%) |
| pplacer | 3 | 26.4 | 5.0% |
| pplacer-scampp | 3 | 188.1 | 6.4% |
| epa-ng-scampp | 3 | 102.4 | 6.9% |
| bscampp | 3 | 8.35 | 5.3% |
| epa-ng | 3 | 10.5 | 5.5% |
| apples2 | 3 | 7.27 | 58.5% (flag: >10%) |
| dipper-m | 3 | 3.14 | 0.6% |
| dipper-m-gpu | 3 | 4.27 | 9.0% |
| dipper-r | 3 | 4.45 | 11.8% (flag: >10%) |

Every arm over ~10 s of non-SEPP placement stays under 7%. The flagged arms either run for a few seconds (apples2, dipper-r: startup and Lustre reads dominate; one apples2 run took twice the others) or are SEPP, which spawns many subprocesses and writes temp files. Flags are reported, not gated.

### Bfull\* go/no-go

Projection as pre-registered: cost linear in backbone size (x32.7 = 327,078 / 10,000) and, for wall time, linear in query count fitted from gg2 V4-150 (1k) and emp-10k (10k) and extrapolated to emp-100k. Peak memory projects the largest B10k\* cell (gg2 full-length) and, separately, the largest V4-150 cell. Limit: 0.8 x 514,000 MB = 402 GB.

| arm | projected emp-100k place wall (h) | projected peak, any region (GB) | projected peak, V4-150 (GB) | verdict (pre-registered rule) |
|---|---|---|---|---|
| sepp-prod | 18.0 | 3,269 | 390 | go on time; memory over 0.8 x node |
| sepp-pinned | 14.1 | 3,297 | 388 | go on time; memory over 0.8 x node |
| pplacer | 21.5 | 473 | 113 | capped attempt only (by design) |
| pplacer-scampp | 98.2 | 1,323 | 816 | no-go on time; memory over 0.8 x node |
| epa-ng-scampp | 52.0 | 1,205 | 1,080 | go, capped on time; memory over 0.8 x node |
| bscampp | 3.0 | 669 | 309 | go on time; memory over 0.8 x node |
| epa-ng | 1.3 | 505 | 201 | B10k\* only (by design) |
| apples2 | 3.5 | 28 | 28 | go |
| dipper-m | 4.1 | 8 | 8 | go |
| dipper-m-gpu | 0.7 | 20 | 20 | go; GPU memory: see probe |
| dipper-r | 6.5 | 15 | 15 | B10k\* only (by design) |

Reading the table:

- The linear-in-backbone memory assumption is known to overstate the decomposition arms: SEPP, pplacer-SCAMPP, EPA-ng-SCAMPP and BSCAMPP place into fixed-size subsets, and production SEPP already runs on the full GG2 reference. Their memory "over 0.8 x node" verdicts mean "may OOM", not "drop". `sepp-prod` is the reference arm and is never dropped.
- Recommendation for T22 (user decision): keep every arm at Bfull\* except the by-design B10k\*-only ones (`epa-ng`, `dipper-r`) and the capped pplacer attempt; an OOM or timeout becomes a DNF row (failures are data). Send full-length cells of the memory-flagged arms to `highmem` (b-cn-25/26, 1.5 TB, 14 days, `AllowAccounts=ALL`; CPU model unchecked, so their wall times carry a hardware caveat).
- `pplacer-scampp` projects 98 h for emp-100k, above the 71 h harness timeout: expect a DNF for that cell unless T22 caps it.
- `epa-ng-scampp` projects 52 h: "go, capped" (a timeout is still a DNF row).

### Findings

- apples2 V4-150 `convert` is DNF on both datasets: `guppy tog` aborts with "Internal edge became -0.0021..." because apples2 writes a `distal_length` longer than the edge. `evaluate` reads the jplace, so accuracy is unaffected; only apples2's V4-150 end-to-end time is missing.
- pplacer `cpu_s` is undercounted (emp-10k: 20 s CPU over 239 s wall): the harness takes CPU from `wait4` rusage, which misses `-j` workers that pplacer does not reap. Wall time and cgroup memory are unaffected.
- Short prep and convert steps fall back to `psutil` when the tool never lifts the cgroup's `memory.peak` above the harness's own baseline (e.g. 13 MB `hmmalign`); all place rows are `cgroup`.
- Zero-length internal edges: B10k\* 0, Bfull\* 1. CPU `dipper-m` passed the RF check at B10k\*; at Bfull\* it may fail (no `--force-binary` in 0.1.5). GPU image `--help` confirms `--force-binary` and no `-T`.
- After `envs/phylobench.yaml` changed, Snakemake scheduled the whole GG2/EMP data build again (software-env provenance). Fixed by `snakemake --cleanup-metadata` on every existing file under `data/` plus `results/emp-*/regions/V4-150/queries.fasta`; T22 must do the same after any change to that env.
- The rejected first stage-2 submissions (`MaxMemPerCPU`) never ran and left no rows.

### Bfull\* DIPPER GPU probe

Pending.

### Actual node-hours

Pending (from `sacct`).
