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

## Results (to be filled after the pilot)

### Per-arm runtime and memory at B10k\*

| arm | dataset | region | wall_mean (s) | peak_mem_max (MB) |
|---|---|---|---|---|
| sepp-prod | gg2 | full-length |  |  |
| sepp-prod | gg2 | V4-150 |  |  |
| sepp-prod | emp-10k | V4-150 |  |  |
| sepp-pinned | gg2 | full-length |  |  |
| sepp-pinned | gg2 | V4-150 |  |  |
| sepp-pinned | emp-10k | V4-150 |  |  |
| pplacer | gg2 | full-length |  |  |
| pplacer | gg2 | V4-150 |  |  |
| pplacer | emp-10k | V4-150 |  |  |
| pplacer-scampp | gg2 | full-length |  |  |
| pplacer-scampp | gg2 | V4-150 |  |  |
| pplacer-scampp | emp-10k | V4-150 |  |  |
| epa-ng-scampp | gg2 | full-length |  |  |
| epa-ng-scampp | gg2 | V4-150 |  |  |
| epa-ng-scampp | emp-10k | V4-150 |  |  |
| bscampp | gg2 | full-length |  |  |
| bscampp | gg2 | V4-150 |  |  |
| bscampp | emp-10k | V4-150 |  |  |
| epa-ng | gg2 | full-length |  |  |
| epa-ng | gg2 | V4-150 |  |  |
| epa-ng | emp-10k | V4-150 |  |  |
| apples2 | gg2 | full-length |  |  |
| apples2 | gg2 | V4-150 |  |  |
| apples2 | emp-10k | V4-150 |  |  |
| dipper-m | gg2 | full-length |  |  |
| dipper-m | gg2 | V4-150 |  |  |
| dipper-m | emp-10k | V4-150 |  |  |
| dipper-m-gpu | gg2 | full-length |  |  |
| dipper-m-gpu | gg2 | V4-150 |  |  |
| dipper-m-gpu | emp-10k | V4-150 |  |  |
| dipper-r | gg2 | full-length |  |  |
| dipper-r | gg2 | V4-150 |  |  |
| dipper-r | emp-10k | V4-150 |  |  |

### Timing CV (3 runs, gg2 V4-150)

| arm | n_ok | wall_mean (s) | wall_cv |
|---|---|---|---|
| sepp-prod |  |  |  |
| sepp-pinned |  |  |  |
| pplacer |  |  |  |
| pplacer-scampp |  |  |  |
| epa-ng-scampp |  |  |  |
| bscampp |  |  |  |
| epa-ng |  |  |  |
| apples2 |  |  |  |
| dipper-m |  |  |  |
| dipper-m-gpu |  |  |  |
| dipper-r |  |  |  |

### Bfull\* go/no-go

| arm | projected emp-100k place wall (h) | projected peak (GB) | verdict |
|---|---|---|---|
| sepp-prod |  |  |  |
| sepp-pinned |  |  |  |
| pplacer |  |  |  |
| pplacer-scampp |  |  |  |
| epa-ng-scampp |  |  |  |
| bscampp |  |  |  |
| epa-ng |  |  |  |
| apples2 |  |  |  |
| dipper-m |  |  |  |
| dipper-m-gpu |  |  |  |
| dipper-r |  |  |  |

### Actual node-hours

| stage | node-hours |
|---|---|
