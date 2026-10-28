# Tickets

One ticket per PR. Move completed tickets to **## Done**. Ticket `T{N}` ↔ GitHub issue `#{N}` in the `phylobench` Project. Design rationale lives in [`design.md`](design.md) (v1) — read the cited sections before planning.

Numbers are never reused or renumbered. Merged or deferred tickets keep their number under **## Merged / Deferred** (issue labeled `deferred` or closed as merged).

Each ticket has: **Week**, **Description**, **Deliverables**, **Acceptance** (what QA checks), **Dependencies**, **Agent notes**.

**Wrapper DoD** (T10, T12, T14, T15, T26):
- [ ] `envs/<tool>.yaml` pinned per `design.md` §1.1 + lockfile (linux-64)
- [ ] `workflow/rules/<tool>.smk` + `config/arms/<arm>.yaml` — **add files only**; never edit shared config (T8 globs them; contract in [`docs/harness.md`](docs/harness.md) §7)
- [ ] Every tool call wrapped by `phylobench-measure` (T4); output normalized to jplace (or T25 adapter) and scored by T5
- [ ] Tool completes on the fixture with T6 ref-prep outputs, in the CI smoke matrix
- [ ] `docs/tools/<tool>.md`: exact command line, version, threads, license of executed engine, quirks

**Local vs CI gates:** agents run on macOS — locally run only `pytest`, `ruff check .`, `ruff format --check .`, `snakemake --lint`. The smoke test (real tool binaries) runs only in CI on `ubuntu-latest`; check it via `gh pr checks`. Unit tests must not require tool binaries.

---

## Open

### Week 2 — Ref-prep, alignment, smoke test, datasets

---

### Week 3 — Real data and pilot

---

### Week 5 — Report

---

### Post-T24 — Non-V4 regions and BSCAMPP follow-up

### T28: Column-window regions (non-V4 support) — closes #59
- **Week:** post-T24
- **Description:** In T22 the six non-V4 regions keep only 46 to 921 of 1,000 queries (V1-V2 153, V1-V3 144, V7-V9 46, against 879 to 921 for V3-V4, V4-150, V4-V5, V5-V6, V6-V8): 711 of the 1,000 sources are MJ0* UMI amplicons without the 27F/1492R primer sites. Primer-trimming of those reads is a hypothesis, UNVERIFIED until the S0 check of `exclusion_counts.tsv` on the cluster. T28 adds a column-window variant of each amplicon region (`V4-150c`, `V3-V4c`, `V1-V2c`, `V1-V3c`, `V4-V5c`, `V5-V6c`, `V6-V8c`, `V7-V9c`) that cuts each query by alignment columns from its shipped aligned row, so a query without primers but covering the window is kept. Labelled the **column-window proxy**: real amplicons start and stop at primers. The primer-matched regions and every T22 output stay untouched and are never pooled with the `c` regions.
- **Window rule (pre-registered):** per region, over the queries' shipped aligned rows that have both primer hits (`regions.amplicon_span`, so V4-150c is the same 150 nt after 515F), window = [median start column, median end column) (`median_low`, integer columns) of the primer-delimited amplicon projected onto columns; window length = median residue count inside the window over those rows.
- **Tolerance (pre-registered, constant across regions):** a query is kept when its non-gap residues inside the window are >= 0.9 x window length (`COVERAGE_TOLERANCE`); the rest are excluded and counted in `exclusion_counts.tsv` as `window coverage below tolerance`.
- **Deliverables:**
  - [x] `regions.column_window`, `regions.cut_column_window`, `COLUMN_REGIONS`; `extract_regions.py --aligned`; region wildcard and schema enum extended (`align.smk`, `config/schema.yaml`)
  - [x] `config/t28.yaml` + `t28-b10k.yaml` (S1) + `t28-bfull.yaml` (S2): sepp-prod, sepp-pinned, bscampp (hmmalign), 8 `c` regions, gg2 only; 24 cells per backbone, 48 in all; T22 rows are reused, nothing is rerun
  - [x] `docs/gg2.md` section, `docs/runbook-t28.md` (S0 exclusion-count check, S1, S2, collect, gates, checksum + copy)
  - [x] Unit tests: window, tolerance edges, exclusion counts, primerless-but-covering query kept, V4-150c length, originals unchanged; overlay cell counts and dry-run DAG
  - [ ] Cluster runs S0, S1, S2 (user runs every command)
- **Acceptance:** `pytest`, `ruff`, `snakemake --lint` and the smoke dry run pass; each T28 dry run plans only `c` cells; `check_stage.py` passes per stage; exclusion counts per `c` region are reported next to the primer-matched ones.
- **Pre-registered reporting criteria (defaults, user may change before S1):**
  - The 0.25-edge mean gate (design section 6) is reported per region and backbone and is stated up front as undecidable at Bfull\* (about 11k to 141k paired queries would be needed for a +-0.25 half-width, against under 1,000 per region).
  - Co-primary: fraction of queries within 3 edges of the truth, bscampp non-inferior to sepp-prod by delta = 2 percentage points (paired per-query bootstrap, upper CI). Labelled **exploratory, chosen after T22**: k = 3 was picked after seeing T22 data, and the sources are reused.
  - Per-region verdicts with CIs; "undecidable" is an allowed outcome; no all-regions rule. Every result is exploratory while the 1,000 sources are reused.
- **Dependencies:** T22, T24
- **Agent notes:** `c` regions apply to gg2 only (EMP has one region). The oracle alignment is not planned for `c` regions. Backbone-independent windows: region queries are dataset-level, so the window comes from the queries' own aligned rows, not a backbone's.

### T29: Accuracy levers at B10k\* — closes #60
- **Week:** post-T24
- **Description:** (a) UPP alignment (bioconda sepp 4.5.6, already pinned in `envs/sepp.yaml`) as aln `upp`, bscampp only (other arms excluded), measured `align` phase, masked to backbone columns (UPP may insert columns); flags and fixed-backbone support verified from the installed package in an S0 probe first; (b) `bscampp-b5000` (`-b 5000 -V 5`) with hmmalign only; b5000 x upp only if both show an effect; optional -V 10 probe. `evaluate.smk` is parameterized by aln (it hard-codes hmmalign). No WITCH (bundled binaries), no oracle at Bfull\*. An alignment lever, not production parity.
- **Deliverables / Acceptance:** to be scoped in plan review (UPP rule, masking unit test, e2e cost in `summaries`/`tables`, `docs/tools/upp.md`, `LICENSES.md` UPP GPL-3, CI smoke only if UPP runs on the tiny fixture in under a minute).
- **Dependencies:** T28 (regions), T22
- **Agent notes:** B10k\* only until a probe sizes UPP at 320k leaves.

### T30: Bfull\* confirmation (conditional on T28/T29) — closes #61
- **Week:** post-T24
- **Description:** Run only if T28/T29 leave a candidate. Confirmatory only on a fresh held-out draw (new backbones; sepp-prod and sepp-pinned rerun) with a criterion fixed before the run; otherwise results stay descriptive. Challenger place/align at Bfull\* routed to `highmem` too, with `cpu_model` recorded and compared only within one CPU model (T22: SEPP rows ran on Xeon 6242, challengers on EPYC 7302).
- **Deliverables / Acceptance:** to be scoped after T28/T29.
- **Dependencies:** T28, T29

---

## Merged / Deferred

- **T9: pplacer wrapper** — merged into T10.
- **T11: pplacer-SCAMPP wrapper** — merged into T12.
- **T13: BSCAMPP wrapper** — merged into T12.
- **T17: RNAsim dataset** — deferred → optional *T17-lite* (one 10k panel, accuracy only) in W4 if on schedule; first cut otherwise. T16 simulated cross-check always runs.
- **T20: SARS-CoV-2 dataset** — deferred out of this cycle (user decision); high-similarity bins are filled on 16S (design §4).

---

## Done

<!-- Move tickets here once merged. Keep the T{N} prefix. -->

### T24: Report and recommendation — closes #24
- [x] `results.md`: methods, verdict per (region, backbone) with effect sizes and CIs, per-bin descriptive tables (gg2 full-length and V4-150), cost tables (gg2, EMP 10k/100k), EMP concordance (labelled concordance), flagged licenses table from `LICENSES.md`; every claim cites a generated table
- [x] Tables are generated by `workflow/scripts/report_tables.py` from the checksum-verified `results/` tables (`src/phylobench/report.py`, pure and tested) and pasted between `<!-- table:key -->` markers; the header carries the regeneration command and the three sha256s
- [x] `tests/test_report.py`: unit tests on the sample tables plus a staleness test (checksums and table blocks vs `results/`), skipped with a reason when the real tables are absent (CI); run locally with `PHYLOBENCH_RESULTS`
- [x] Caveats: alpha=1.0, pruned Bfull\*, uDance/APPLES-2 truth bias and the sim10k cross-check (B10k\* only), DIPPER joint placement and zero-length contraction, hmmalign arm (oracle gap), short window, `mem_source`, hardware split, DNFs, B10k\*-only arms, small regions
- [x] Recommendation from accuracy, runtime and memory only: no arm passes the gate at Bfull\* (effect sizes, `sepp-pinned` as noise floor, no equivalence claim); B10k\* gate and cost ranking; GPU beside the CPU ranking; licenses flagged, not a factor
- [x] User sign-off (PR review, 2026-10-28)

### T23: Benchmark plots — closes #23
- [x] `tests/fixtures/sample_results/{runs,accuracy}.tsv` (about 240 KB, synthetic, schema-exact; gg2 + sim10k, 11 arms, apples2 convert DNFs, all-unplaced pplacer at Bfull\*, a timeout DNF, psutil and Xeon rows, GPU memory)
- [x] 10 PNGs: accuracy (edge distance, exact fraction, delta error) by similarity bin, backbone x region; verdict (paired difference vs sepp-prod, bootstrap 95% CI, pass margin 0.25, `unplaced`/`no-reference` states, sim10k cross-check); oracle-vs-hmmalign gap; runtime and peak memory vs similarity, backbone size and query count (log axes, DNF `x`, psutil hollow); end-to-end vs placement-only; GPU memory
- [x] `src/phylobench/summaries.py` (pure, tested; reused by T24) and `workflow/scripts/plot_{accuracy,cost}.py` driven only by `runs.tsv`/`accuracy.tsv`
- [x] Cost vs similarity uses each cell's median query similarity (captioned): runtime is measured per cell, not per bin
- [x] `snakemake plots -c 1` regenerates all figures from the sample TSVs in CI; real tables are plotted with the scripts directly (README)

### T22: Full benchmark execution — closes #22
- [x] `config/full.yaml` (304 cells; `exclude` key drops B10k\*-only arms, oracle and sim10k at Bfull\*, one capped pplacer Bfull\* cell) staged as S1 `full-b10k.yaml` (215 cells, `profiles/slurm-pilot`) and S2 `full-bfull.yaml` (89 cells, `profiles/slurm`); S3 collect with `full.yaml` alone; `check_stage.py` gates (#52)
- [x] Routing from measured memory: SEPP place rules on `highmem` (Intel Xeon Gold 6242, hardware caveat), every other place rule on whole `long` nodes (AMD EPYC 7302) (#54)
- [x] DIPPER adapter contracts zero-length backbone edges (CPU dipper-m scored at Bfull\*) (#53); pplacer `-nan`/`-inf` placements scored as unplaced (#55)
- [x] sim10k panel: 11 arms x {full-length, V4-150}, hmmalign, B10k\* only; T17-lite skipped (user decision)
- [x] Acceptance (`check_stage.py config/full.yaml`, 2026-10-28): coverage, mem_source, cpu_ratio, threads, gpu PASS; 839 `runs.tsv` rows over 304 cells; no OOM or timeout DNFs
- [x] DNFs (all `convert`, data not failures): 29 apples2 `guppy tog` (allow-listed); 1 pplacer at gg2/Bfull\*/V4-150 (`guppy tog` cannot lex `-nan`). That pplacer cell placed ok in 2.6 h / 74 GB but its likelihoods underflowed: 0 of 879 queries scored
- [x] Node-hours: P0 0.44, S1 8.13, S2 51.04 (highmem 22.41, long 28.03, gpu 0.51, short 0.09), S3 local; total about 60 against a 580 ceiling
- [x] Archive at git `8fb2ab8`: `/ddn_scratch/lxxu/phylobench-archive/t22-final.tgz` (958 MB, 3,252 entries, sha256 `234bbbba27ab0b3d609e6690d24a95a4eb0acc31985ce7686a7ec4016128cb37`); tables in `t22-final/` with `SHA256SUMS`, copied to laptop `results/` and verified

### T21: HPC profile and pilot — closes #21
- [x] `profiles/slurm/config.yaml`: one rule per sbatch (cgroup `memory.peak` primary), retries 0, 16 cpus/task, `--exclusive` place rules on `long` (72 h), shared prep/convert, DIPPER GPU on `gpu` / `b-gpu-21` (RTX 3090, `gpu: 1`, 196608 MB under `MaxMemPerCPU`); `profiles/slurm-pilot` moves `long` rules to `short` (higher PriorityTier); executor plugins 2.8.0 / jobstep 0.6.1 pinned; harness `timeout_s` 71 h
- [x] Pilot at B10k\*: all 33 cells (11 arms x gg2 full-length, gg2 V4-150, emp-10k V4-150) `ok`, 3 repeats on gg2 V4-150 for timing CV, plus a Bfull\* `dipper-m-gpu` probe (GPU 780 MB); about 2.1 node-hours
- [x] `docs/pilot.md`: per-arm runtime/memory, timing CV, pre-registered Bfull\* go/no-go with T22 recommendation (keep all arms, `highmem` for memory-flagged full-length cells), findings (apples2 `guppy tog` DNF, pplacer `cpu_s` undercount, env-provenance metadata cleanup)

### T27: Day-0 cluster capability check — closes #27
- [x] `docs/cluster.md`: Barnacle nodes (GPU node `b-gpu-23`: 2× Xeon Silver 4110, 16 CPUs, 93 GiB, kernel 6.12 EL10), account `knightlab`, partitions and wall limits, scratch `/ddn_scratch/$USER` (no quota)
- [x] cgroup v2 `memory.peak` readable from a job step → `mem_source=cgroup` primary
- [x] Apptainer 1.5.0, driver 580.126.20 / CUDA 13.0; RTX 3090 on `b-gpu-21`/`-22`, Tesla V100 32 GB on `b-gpu-23`/`-24` (probed); pinned DIPPER image runs a fixture `--add` with `--nv` on both cards (no `-T` in the GPU build)
- [x] `snakemake-executor-plugin-slurm` 2.8.0 (+ jobstep 0.6.1) resolves against Snakemake 9.27.0; pinning is T21
- [x] Memory sections of `deployer.md` / `diagnostics.md`

### T18: GG2 16S dataset — closes #18
- [x] Checksummed download via `workflow/scripts/fetch.sh` (shared with T19); both `.qza` sha256 pinned in `workflow/rules/gg2.smk`; `.qza` unzipped with stdlib `zipfile` (`unzip_qza.py`)
- [x] Bfull\* = shipped minus held-out minus exact duplicate rows (0 on GG2 2022.10); B10k\* = seeded uniform subsample; `check_backbone` enforces bifurcation, no held-out leaf, no duplicates (`src/phylobench/gg2.py`)
- [x] 1,000 held-out queries × 9 regions; ≥99.9 seeding by replacement; `data/gg2/MANIFEST.tsv` (`kind backbone region bin path n sha256 flagged`), bins < 100 flagged; L = 1,598, R = 331,269 (`docs/gg2.md`)
- [x] `accuracy.true_positions` accepts extra non-backbone reference leaves (primer exclusions)
- [ ] Full build on HPC (user-confirmed): `snakemake --use-conda -c 16 data/gg2/MANIFEST.tsv`

### T19: EMP V4 ASV query sets — closes #19
- [x] Checksummed BIOM download (`emp_download` via T18's `fetch.sh`; sha256 pinned on first download, none published); all 197,711 observation IDs verified unique 150-bp ACGT
- [x] Nested `emp-10k` / `emp-100k` query sets (seed 0), region `V4-150`, `alns: [hmmalign]`; `sim` per ASV vs each backbone; `data/emp/MANIFEST.tsv` (T18 columns + `dataset`)
- [x] Concordance = edge distance to the `sepp-prod` placement on the same backbone (`phylobench.concordance`, `concordance` rule) into `results/concordance.tsv`, `label=concordance`; `accuracy.tsv` carries accuracy only
- [x] EMP links GG2's backbone, ref-prep and HMM (one per backbone); `config/emp.yaml`; `docs/datasets/emp.md`; `h5py` added to `envs/phylobench.yaml`

### T12: bscampp family + EPA-ng control — closes #12
- [x] Wrapper DoD for `pplacer-scampp`, `epa-ng-scampp`, `bscampp`, `epa-ng`: `workflow/rules/{bscampp,epa-ng}.smk`, `config/arms/*.yaml`, `docs/tools/{bscampp,epa-ng}.md`
- [x] External pinned `pplacer`/`epa-ng`, per-env config from `envs/bscampp.post-deploy.sh`; rules check `place.log` for the env config paths and fail if `~/.bscampp` is initialized
- [x] `bscampp`: `--threads {threads} --cpus-per-job 2`, subtree 2000, votes 5 explicit
- [x] `epa-ng` control is B10k\*-only by design; not expressible in arm config, so T22 drops it at Bfull\*

### T10: pplacer + sepp-pinned wrappers — closes #10
- [x] `pplacer` arm at Wrapper DoD (`-j {threads}`, shipped `raxml-info.txt`); harness `--timeout 86400` on every pplacer cell (deviation from design.md §5: global, not Bfull\*-only)
- [x] `sepp-pinned` arm (required: QIIME 2 2026.7 ships pplacer alpha19, pin is alpha23): `workflow/rules/sepp-pinned.smk`, `config/arms/sepp-pinned.yaml` — `run_sepp.py -A 1000 -P 5000 -x {threads}`, cell-local `-d`/`-p`, unaligned queries, phases `align+place` + shared `convert`; `smoke_min_exact` to be set once CI reports it
- [x] Capped Bfull\* overlay `config/pplacer-bfull-cap.yaml` (placeholder dataset `gg2-v4-1k`, produced by T16/T21 data prep)
- [x] SEPP's jplace v1 normalized to v3 by T26's `sepp_jplace.py` (`phylobench.jplace.from_sepp`) after an ok run
- [x] Provenance: identical pplacer url#md5 lock line across pplacer/sepp/bscampp locks, enforced by `envs.yml` lock-guard; `env_lock_sha256` on every row
- [x] `docs/tools/pplacer.md`, `docs/tools/sepp-pinned.md`; `LICENSES.md` unchanged

### T26: sepp-prod wrapper (production baseline) — closes #26
- [x] Run SEPP standalone with production's command and versions (QIIME 2 2026.7's sepp/pplacer/hmmer) on the GG2 reference, outside QIIME 2 — the yardstick for every verdict.
- [x] Wrapper DoD (env = `envs/sepp-prod.yaml`: sepp 4.5.6, pplacer 1.1.alpha19, hmmer 3.4, openjdk 25.0.2 as shipped in QIIME 2 2026.7, own lockfile; no QIIME 2)
- [x] `run-sepp.sh <q> <name> -x {threads} -A 1000 -P 5000 -a -t -r` directly, per-job working dir and temp dir; single phase `align+place+convert`
- [x] SEPP's jplace (v1, `[N]` edges) normalized to v3 and scored by the shared `evaluate` rule
- Scope change 2026-10-19 (user): no q2-fragment-insertion, no QIIME 2 env, no q2-vs-direct equality check.
- Files: `workflow/rules/sepp-prod.smk`, `config/arms/sepp-prod.yaml`, `envs/sepp-prod.yaml`, `phylobench.jplace.from_sepp`, `docs/tools/sepp.md`

### T15: DIPPER wrapper — closes #15
- [x] Arms `dipper-m` (CPU, `-i m -d 2 -T`), `dipper-r` (raw `-i r -T`, B10k\* only via run config), `dipper-m-gpu` (Apptainer `--nv`, harness `--gpu`, image pinned by digest; manual tier, verified on the cluster in T27)
- [x] `workflow/rules/dipper.smk` → T25 adapter (`dipper_to_jplace.py`, fail-loud RF check) → T5 evaluate; CPU arms green in the CI smoke matrix
- [x] `docs/tools/dipper.md`, `LICENSES.md`; CI smoke glob skips `*-gpu` arms

### T8: Tiny end-to-end smoke test — closes #8
- [x] `config/smoke.yaml`, `config/arms/pplacer.yaml`, `workflow/rules/{pplacer,evaluate,smoke}.smk`; `phylobench.tables` + `check_smoke.py`
- [x] `results/runs.tsv` (design §6 columns + `run_id`, `dnf_reason`, `env_lock_sha256`; non-ok → `status=DNF`) and `results/accuracy.tsv`
- [x] CI `arms` → `smoke-arm` matrix (one job per `config/arms/*.yaml`) → `smoke` aggregate; asserts backbone-width aligned queries, forced-DNF row, ok + provenance on every other row, no conda pin-file fallback
- [x] Arm contract for wrappers: `docs/harness.md` §7; `docs/tools/pplacer.md`; guppy layout `{cell}/placements.tog.*`

### T14: APPLES-2 wrapper — closes #14
- [x] `workflow/rules/apples.smk` (`apples2_place`, `run_apples.py -D` on the T6 `me-tree.nwk`, no per-job FastTree), `config/arms/apples2.yaml`
- [x] `docs/tools/apples.md`: output numbering, uDance-bias caveat, known unplaceable-query issue (T22)

### T7: Regions and query alignment — merged in #41 (2026-10-16)
`phylobench.regions` (IUPAC primer matching for V4-150, V3-V4, V1-V2, V1-V3, V4-V5, V5-V6, V6-V8, V7-V9 with exclusion counts) and `workflow/rules/align.smk` (backbone `hmmbuild` + `hmmalign`, masked to backbone columns, timed `align` phase, oracle passthrough for B10k\*). Notes in `docs/tools/hmmer.md`.

### T16: Similarity, held-out construction, and simulated cross-check — merged in #39 (2026-10-16)
`phylobench.similarity` (windowed exact-integer sim, bins, per-bin report), `phylobench.heldout` (monophyletic-taxon clade/singleton selection, tree/alignment pruning), `workflow/scripts/make_sim_panel.py` (10k-leaf simulated panel to `data/sim10k`). ≥99.9 seeding moved to T18 (needs Bfull* similarity). 10k × 10k × 1500 similarity ≈ 5 s.

### T25: Tree <-> placement conversion (guppy tog + DIPPER adapter)
- [x] `workflow/rules/convert.smk`: `guppy tog` rule, harness-wrapped (phase `convert`)
- [x] `phylobench.dipper_adapter`: XOR edge signatures, one postorder pass, topology check fails loudly; reuses T5 `true_positions`
- [x] `guppy-tog` CI workflow (`convert.yml`, not required); adapter perf test on 100k queries / 330k leaves < 60 s
- [x] `docs/tools/dipper.md`, `docs/tools/guppy.md`; design.md 1.3 rooting-agnostic

### T2: License and install audit — merged in #33 (2026-10-15)
Pinned envs + linux-64 lockfiles for 8 tools and QIIME 2 2026.7, `envs.yml` CI (lock-guard + 9 installs, not required), `LICENSES.md` (flags only), `docs/audit.md`. pplacer pinned at 1.1.alpha23; QIIME 2 2026.7 ships alpha19 → `sepp-pinned` required. BSCAMPP uses an immutable per-env config.

### T3: Tiny fixture dataset — merged in #31 (2026-10-15)
Platform-deterministic simulator (`phylobench.simulate`), 100-leaf backbone + 10 held-out queries in `tests/fixtures/tiny/`, byte-identical on macOS and Linux.

### T4: Measurement harness — merged in #32, test fix #35 (2026-10-15)
`phylobench-measure` / `python -m phylobench.measure`; DNF and invocation contract in `docs/harness.md`. Linux `ru_maxrss` is inherited across fork/exec — never call `run()` in-process from a large parent.

### T5: Accuracy metrics library — merged in #34 (2026-10-15)
`phylobench.jplace` + `phylobench.accuracy`: bipartition-keyed truth, polytomy-aware edge distance, delta error, public `edge_distance` core for T19, hashing helpers for T25. 331k-leaf scoring ≈ 17 s.

### T1: Repo skeleton and CI — merged in #29 (2026-10-15)
Skeleton, `config/schema.yaml`, `envs/phylobench.yaml`, CI job `lint-test` (now a required check on `main`). Pin change: pandas 2.3.3 (bioconda snakemake 9.27.0 requires pandas <3).

### T6: Reference prep (pruning + model conversion) — merged 2026-10-16
`phylobench.models` (RAxML-8 info to `GTR{..}+FU{..}+G4m{1.0}`, strict `ac ag at cg ct gt` label check, alpha as shipped), `refprep.smk` (model writer + harness-wrapped FastTreeMP `-nome -noml` ME tree, phase `ref-prep`), `envs/raxml-ng.yaml` 2.0.3 + lock. `envs.yml` T6 steps: raxml-ng `--check`, epa-ng and pplacer on the fixture, FastTreeMP run. No pruning (T16/T18).
