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

### Week 1 — Foundations

### T27: Day-0 cluster capability check
**Week:** W1 (day 0) · **Owner:** human + Deployer (needs cluster access)
**Description:** Confirm the HPC environment supports the measurement and GPU plan before code depends on it.

**Deliverables:**
- [ ] `docs/cluster.md`: node CPU model/cores/RAM, kernel version, SLURM account/partition, scratch quota, wall limit
- [ ] Whether a SLURM job can read its own cgroup v2 `memory.peak` (`srun cat /sys/fs/cgroup/$(...)/memory.peak`)
- [ ] `apptainer --version`; NVIDIA driver/CUDA version on the RTX 3090 node; pinned DIPPER image runs `--help` under Apptainer with `--nv`
- [ ] Snakemake SLURM executor plugin availability
- [ ] Fill `diagnostics.md` / `deployer.md` Memory with hostnames and access pattern

**Acceptance:** each item answered yes/no with command output pasted (no secrets).

**Dependencies:** none

**Agent notes:** No job submission beyond trivial `srun` probes. Decides T4's primary memory source.

---

### Week 2 — Ref-prep, alignment, smoke test, datasets

---

### T18: GG2 16S dataset
**Week:** W2
**Description:** Download the GG2 2022.10 SEPP reference + taxonomy and build Bfull\*, B10k\*, and the held-out query sets for all nine regions.

**Deliverables:**
- [ ] Scripted checksummed download of `2022.10.backbone.sepp-reference.qza` and `2022.10.backbone.tax.qza`; unzip `.qza` without needing QIIME 2
- [ ] Bfull\* (shipped minus held-out, via T16) and B10k\* (uniform subsample of Bfull\*); **pruned-reference checksums recorded** in `data/gg2/MANIFEST.tsv`
- [ ] 1,000 held-out source sequences × {full-length, V4-150, V3–V4, V1–V2, V1–V3, V4–V5, V5–V6, V6–V8, V7–V9}; per-size bin counts; ≥99.9 seeding applied

**Acceptance:** manifest lists counts per backbone × region × bin and checksums; no held-out leaf in either backbone (test); bins < 100 flagged.

**Dependencies:** T7, T16

**Agent notes:** design §3–§4. `data/` never committed. ≥99.9 seeding (design §4 step 2) built here from similarity() + heldout.choose_groups; moved from T16. Similarity memory ≈16·R·L bytes (R reference leaves, L alignment columns); check against GG2 alignment width.

---

### Weeks 2–3 — Tool wrappers (parallel; add files only)

### T10: pplacer + sepp-pinned wrappers
**Week:** W2–W3
**Description:** pplacer baseline arm on the shared pinned binary, plus the conditional `sepp-pinned` arm.

**Deliverables:**
- [ ] Wrapper DoD for `pplacer` (`-j {threads}`, shipped `raxml-info.txt`)
- [ ] `sepp-pinned` only if `docs/audit.md` (T2) says the QIIME 2 versions differ from pins; otherwise document "not needed" and stop
- [ ] Capped Bfull\* pplacer config: 1k V4 queries, harness `--timeout 86400`

**Acceptance:** Wrapper DoD; provenance shows the same pplacer binary hash as other pinned arms.

**Dependencies:** T8

**Agent notes:** alpha23 pin (alpha20 deadlocks with `-j`; alpha22→23 build-only).

---

### Week 3 — Real data and pilot

### T19: EMP V4 ASV query sets
**Week:** W3
**Description:** Real EMP Deblur V4 ASVs for runtime/memory scaling and concordance with sepp-prod.

**Deliverables:**
- [ ] Scripted checksummed download of `ftp.microbio.me/emp/release1/otu_tables/deblur/emp_deblur_150bp.release1.biom`; verify observation IDs are sequences
- [ ] Query sets of 10k and 100k unique ASVs; per-ASV `sim` vs. each backbone
- [ ] Concordance = edge distance to `sepp-prod` placement on the same backbone, labeled `concordance`

**Acceptance:** manifest; no output column labels concordance as accuracy.

**Dependencies:** T18, T26

**Agent notes:** GG2 ASV positions came from DEPP — not truth.

---

### T21: HPC profile and pilot
**Week:** W3
**Description:** SLURM profile and a pilot run to calibrate resources and decide per-arm go/no-go at Bfull\*.

**Deliverables:**
- [ ] `profiles/slurm/config.yaml`: exclusive nodes, 16 cpus/task, node memory, 72 h cap, GPU partition for DIPPER
- [ ] Pilot at B10k\*: dipper-m (GPU, CPU) and sepp-prod first, then each merged arm
- [ ] `docs/pilot.md`: per-arm runtime/memory at B10k\*, timing CV (3 repeats), go/no-go for Bfull\*

**Acceptance:** pilot `runs.tsv` covers every merged arm at B10k\*; go/no-go table complete.

**Dependencies:** T8, T18, T26, T27

**Agent notes:** Deployer submits only after the user confirms account and node-hours. Other wrappers join as they merge.

---

### T23: Benchmark plots
**Week:** W3 (scripts) → W5 (final)
**Description:** Figure scripts built against sample TSVs, then run on real results.

**Deliverables:**
- [ ] `tests/fixtures/sample_results/{runs,accuracy}.tsv` (small, synthetic, schema-exact)
- [ ] Accuracy (edge distance, exact fraction, delta error) vs. similarity bin, faceted by backbone × region; verdict panel (paired difference vs. sepp-prod with CIs)
- [ ] Runtime and peak memory vs. similarity bin, vs. backbone size, and vs. query count (log–log, DNF markers); end-to-end vs. placement-only; GPU memory panel; oracle-vs-hmmalign gap
- [ ] `workflow/scripts/plot_*.py` driven only by `results/*.tsv`

**Acceptance:** `snakemake plots` regenerates all figures from the sample TSVs in CI.

**Dependencies:** T8 (scripts); T22 (final figures)

**Agent notes:** Never pool bins across backbone sizes.

---

### Week 4 — Full runs

### T22: Full benchmark execution
**Week:** W4
**Description:** Staged runs on HPC: B10k\* accuracy → Bfull\* accuracy → EMP scaling.

**Deliverables:**
- [ ] `config/full.yaml` from T21 go/no-go
- [ ] Complete `results/runs.tsv` + `results/accuracy.tsv`, DNFs recorded
- [ ] Archived results tarball with checksums and env locks
- [ ] Accuracy arms on the T16 simulated 10k panel (data/sim10k) — uDance/APPLES-2 cross-check

**Acceptance:** every planned (backbone, region, arm) cell has a row with status — no silent gaps.

**Dependencies:** T21, T18, T19, T10, T12, T14, T15, T26

**Agent notes:** Deployer owns submission; two-failed-attempt rule per job class. B10k\*-only arms (`pplacer`, `epa-ng`, `dipper-r` per `design.md` §5) cannot be expressed in arm config (`config/schema.yaml` has `additionalProperties: false` and `plan_cells` runs every arm on every backbone), so drop them at Bfull\* with a staged `config/full.yaml` or a separate `--configfile`. Optional T17-lite only if B10k\* is done by W4 day 2.

---

### Week 5 — Report

### T24: Report and recommendation
**Week:** W5
**Description:** Results report with the SEPP verdict per region × backbone, cost tables, a flagged license table, and a recommendation.

**Deliverables:**
- [ ] `results.md`: methods summary, verdict per (region, size) with per-bin descriptive results, runtime/memory tables, license table from `LICENSES.md` (flag only)
- [ ] Caveats: alpha=1.0 model, pruned Bfull\*, uDance/APPLES-2 truth bias + simulated cross-check, DIPPER joint placement, hmmalign alignment arm, short-window ambiguity, `mem_source`
- [ ] Recommendation based on accuracy, runtime, and memory only; GPL/AGPL noted as flags, not factors

**Acceptance:** every claim cites a figure or table; user sign-off.

**Dependencies:** T23, T2

**Agent notes:** design §0, §9. No license analysis.

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

### T12: bscampp family + EPA-ng control — closes #12
- [x] Wrapper DoD for `pplacer-scampp`, `epa-ng-scampp`, `bscampp`, `epa-ng`: `workflow/rules/{bscampp,epa-ng}.smk`, `config/arms/*.yaml`, `docs/tools/{bscampp,epa-ng}.md`
- [x] External pinned `pplacer`/`epa-ng`, per-env config from `envs/bscampp.post-deploy.sh`; rules check `place.log` for the env config paths and fail if `~/.bscampp` is initialized
- [x] `bscampp`: `--threads {threads} --cpus-per-job 2`, subtree 2000, votes 5 explicit
- [x] `epa-ng` control is B10k\*-only by design; not expressible in arm config, so T22 drops it at Bfull\*

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
