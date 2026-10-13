# Tickets

One ticket per PR. Move completed tickets to **## Done**. Design rationale for every ticket lives in [`design.md`](design.md) — read the cited section before planning.

Each ticket has: **Milestone**, **Description**, **Deliverables**, **Acceptance** (what QA checks), **Dependencies**, and **Agent notes** (gotchas the Planner must address).

Shared Definition of Done for tool-wrapper tickets (T9–T15), referred to as **Wrapper DoD**:
- [ ] `envs/<tool>.yaml` pinned to the exact version in `design.md` §1.1 (linux-64)
- [ ] `workflow/rules/<tool>.smk` rule(s) for `place`, wrapped by the T4 measurement harness
- [ ] Output normalized to jplace (or the shared placement TSV) and scored by T5
- [ ] Tool added to the smoke-test matrix (`config/smoke.yaml`) and green in CI
- [ ] Thread / config pinning documented in `docs/tools/<tool>.md` (command line, version, license, known quirks)

---

## Open

### M0 — Foundations (weeks 1–2)

### T1: Repo skeleton and CI
**Milestone:** M0
**Description:** Create the Snakemake/Python project skeleton and a GitHub Actions CI that runs tests and linters on `ubuntu-latest` (x86_64).

**Deliverables:**
- [ ] Layout: `workflow/Snakefile`, `workflow/rules/`, `workflow/scripts/`, `src/phylobench/` (Python package), `envs/`, `config/`, `tests/`, `docs/tools/`
- [ ] `pyproject.toml` (Python 3.11, pytest, ruff, pandas, matplotlib, treeswift) and `envs/phylobench.yaml`
- [ ] `config/schema.yaml` validating datasets, sizes, tools, arms, threads
- [ ] `.github/workflows/ci.yml`: `ruff check .`, `ruff format --check .`, `snakemake --lint`, `pytest`
- [ ] `.gitignore` adds `data/`, `results/`, `.snakemake/`, `__pycache__/`, `*.sif`
- [ ] README Setup section with install/test/lint commands; fill QA.md placeholders if they diverge

**Acceptance:** CI green on a PR with a trivial passing test; `snakemake -n` succeeds on empty config.

**Dependencies:** none

**Agent notes:** KISS — no tool rules yet. Use mamba/micromamba in CI with env caching.

---

### T2: License and install audit
**Milestone:** M0
**Description:** Install every tool from a pinned env on linux-64, confirm every `[UNVERIFIED]` item in `design.md` §1, and record licenses including transitive/bundled ones.

**Deliverables:**
- [ ] `envs/{pplacer,sepp,bscampp,apples,dipper,epa-ng,raxml-ng,fasttree,raxml8,witch,mafft}.yaml`, each with an explicit lockfile
- [ ] `LICENSES.md`: tool → own license → effective license of executed engine(s) → bundled binaries → source URL + commit/version checked
- [ ] `docs/audit.md`: resolution of every `[UNVERIFIED]` in design §1.2–§1.4 (pplacer model-file formats, SEPP FastTree-log support, APPLES thread flag, DIPPER backbone preservation, TBB license)
- [ ] pplacer pin decision (alpha22 vs. alpha19) with rationale
- [ ] Proof BSCAMPP/SCAMPP can be forced to use pinned external `pplacer`/`epa-ng` and a per-job config (no `~/.bscampp` reuse)
- [ ] DIPPER GPU image pulled by digest and run under Apptainer on the target GPU node (or documented as blocked)

**Acceptance:** each env solves and `<tool> --version` (or equivalent) prints the pinned version in CI; `LICENSES.md` has a row for all 8 placement arms plus ref-prep/alignment tools.

**Dependencies:** T1

**Agent notes:** Read LICENSE files at the pinned tag, not `main`. The question is "what executes", not "what the wrapper repo says". Flag any GPL/AGPL in the execution path.

---

### T3: Tiny fixture dataset
**Milestone:** M0
**Description:** Commit a tiny simulated dataset with known true placements for unit and smoke tests.

**Deliverables:**
- [ ] `tests/fixtures/tiny/`: ~100-leaf true tree, true alignment, 10 full-length + 10 fragment queries, `truth.tsv` (query → true edge)
- [ ] `workflow/scripts/make_fixture.py` that regenerates it deterministically (fixed seed)
- [ ] Total size < 1 MB

**Acceptance:** regeneration script reproduces committed files byte-identically.

**Dependencies:** T1

**Agent notes:** Subsample from RNAsim if quick to fetch; otherwise simulate. Keep sequences long enough (~1,000 sites) that tools don't fail on degenerate input.

---

### T4: Measurement harness
**Milestone:** M0
**Description:** A single wrapper that runs a command and records wall time, CPU time, peak memory, GPU memory, and exit status as JSON.

**Deliverables:**
- [ ] `src/phylobench/measure.py` + CLI `phylobench-measure --out run.json -- <cmd>`
- [ ] Peak memory from cgroup v2 `memory.peak` when available (primary), `/usr/bin/time -v` max RSS (secondary); record which source was used
- [ ] Optional GPU poller (`nvidia-smi` every 0.5 s) behind `--gpu`
- [ ] Status classification: `ok | oom | timeout | error`
- [ ] Provenance fields: node, CPU model, threads, tool version, env lock hash, command

**Acceptance:** unit tests covering status classification and JSON schema; a test that a forking command (e.g. `bash -c 'python -c alloc & python -c alloc; wait'`) reports a peak ≥ the sum when cgroup is available.

**Dependencies:** T1

**Agent notes:** design §6. `time -v` underreports forking tools — never let it be the silent primary. Keep it one module; no plugin system.

---

### T5: Accuracy metrics library
**Milestone:** M0
**Description:** Pure functions to read placements and score them against a true tree.

**Deliverables:**
- [ ] `src/phylobench/jplace.py`: minimal jplace reader → (query, edge_num, LWR) top placement
- [ ] `src/phylobench/accuracy.py`: polytomy-aware edge distance (collapse zero-length branches; any edge in the true polytomy set = 0); delta error (FN bipartitions, Mirarab et al. 2012)
- [ ] Unit tests on hand-built trees, incl. polytomy and zero-length cases

**Acceptance:** 100% branch coverage of `accuracy.py`; metrics agree with hand-computed values.

**Dependencies:** T1

**Agent notes:** design §3.2. Must scale to 200k-leaf trees — test that edge distance on a 200k random tree runs in seconds, not minutes.

---

### T6: Reference-prep rules
**Milestone:** M0
**Description:** Snakemake rules producing each tool family's backbone parameter files on a fixed topology, each timed separately.

**Deliverables:**
- [ ] `raxml-ng --evaluate` GTR+G4 → `.bestModel` (EPA-ng family)
- [ ] FastTree GTR+CAT log (pplacer family)
- [ ] RAxML-8 `-f e` info file (SEPP)
- [ ] FastTree `-nome -noml` ME branch lengths (APPLES-2)
- [ ] All wrapped by T4 harness; phase = `ref-prep`

**Acceptance:** rules run on the T3 fixture in CI; outputs accepted by each downstream tool's `--help`-documented input parser.

**Dependencies:** T2, T4

**Agent notes:** design §2 "Model parameters". Topology is never re-estimated per tool.

---

### T7: Query alignment rules
**Milestone:** M0
**Description:** Produce the shared arm-A query alignment, masked to backbone columns.

**Deliverables:**
- [ ] True-alignment passthrough (RNAsim / fixture)
- [ ] WITCH path (16S)
- [ ] Reference-based path for SARS-CoV-2 (`mafft --addfragments --keeplength`)
- [ ] Mask/trim to backbone width; phase = `align`, timed via T4

**Acceptance:** on the fixture, every aligned query has exactly the backbone's column count; WITCH output identical across two runs with the same seed.

**Dependencies:** T2, T4

**Agent notes:** design §2 "Query alignment". Open alternative: reuse SEPP's extended alignment — raise in review, don't silently choose.

---

### T8: Tiny end-to-end smoke test
**Milestone:** M0
**Description:** Wire fixture → ref-prep → align → pplacer → evaluate → `runs.tsv` + `accuracy.tsv`, running in CI.

**Deliverables:**
- [ ] `config/smoke.yaml` with the fixture and `tools: [pplacer]`
- [ ] `workflow/rules/evaluate.smk` producing tidy `results/runs.tsv` and `results/accuracy.tsv` (columns per design §6)
- [ ] CI job `smoke` running `snakemake --use-conda --configfile config/smoke.yaml -c 2`

**Acceptance:** CI smoke job green; `accuracy.tsv` has 20 rows; pplacer full-length edge error on the fixture is ≤ 1 for ≥ 8/10 queries.

**Dependencies:** T3, T4, T5, T6, T7

**Agent notes:** This is the template every wrapper (T9–T15) copies. Keep the pplacer rule minimal; T9 hardens it.

---

### M1 — Tool wrappers (weeks 2–4)

### T9: pplacer wrapper (baseline)
**Milestone:** M1
**Description:** Harden the pplacer rule into the production baseline arm.

**Deliverables:**
- [ ] Wrapper DoD
- [ ] `-j {threads}`; FastTree-log refpkg input; pinned binary path exported for reuse by T10/T11

**Acceptance:** Wrapper DoD; same pplacer binary hash recorded in provenance for T9/T10/T11 runs.

**Dependencies:** T8

**Agent notes:** design §1.1 pin note (alpha20 deadlocks with `-j`).

---

### T10: SEPP wrapper
**Milestone:** M1
**Description:** Run SEPP natively on unaligned queries, in both the lab's production configuration and default configuration.

**Deliverables:**
- [ ] Wrapper DoD
- [ ] Two arms: `sepp-prod` (config from Q2 answer; default q2-fragment-insertion settings) and `sepp-default`
- [ ] Forced to use T9's pinned pplacer
- [ ] Phase recorded as `align+place` (inseparable)

**Acceptance:** Wrapper DoD; both arms produce jplace on the fixture.

**Dependencies:** T8

**Agent notes:** design §2 "SEPP config". License: GPL-3. SEPP rebuilds HMMs every run — that cost stays in.

---

### T11: pplacer-SCAMPP wrapper
**Milestone:** M1
**Description:** Run SCAMPP with pplacer as the placement engine via `bscampp`'s `run_scampp.py`.

**Deliverables:**
- [ ] Wrapper DoD
- [ ] `--placement-method pplacer`, pinned external pplacer (T9), per-job config dir (no `~/.bscampp`)

**Acceptance:** Wrapper DoD; provenance shows external (not bundled) pplacer path.

**Dependencies:** T8

**Agent notes:** MIT wrapper, GPL-3 engine — note in `docs/tools/`.

---

### T12: EPA-ng-SCAMPP wrapper + plain EPA-ng control
**Milestone:** M1
**Description:** Run SCAMPP with EPA-ng, plus plain EPA-ng as a control arm to separate engine from framework effects.

**Deliverables:**
- [ ] Wrapper DoD for both `epa-ng-scampp` and `epa-ng`
- [ ] Pinned external epa-ng 0.3.8; raxml-ng `.bestModel` input from T6; per-job config

**Acceptance:** Wrapper DoD for both arms; provenance shows external epa-ng.

**Dependencies:** T8

**Agent notes:** EPA-ng is **AGPL-3**. Control arm subject to Q7 (default: include).

---

### T13: BSCAMPP wrapper
**Milestone:** M1
**Description:** Run BSCAMPP (EPA-ng base) with documented defaults.

**Deliverables:**
- [ ] Wrapper DoD
- [ ] `--threads {threads} --cpus-per-job 2`; subtree size 2000, votes 5 (defaults, recorded explicitly); pinned external epa-ng; per-job config

**Acceptance:** Wrapper DoD; concurrent smoke runs don't clobber each other's config.

**Dependencies:** T8

**Agent notes:** design §1.1 BSCAMPP bundled-binary risk.

---

### T14: APPLES-2 wrapper
**Milestone:** M1
**Description:** Run APPLES-2 on aligned queries against an ME-branch-length backbone.

**Deliverables:**
- [ ] Wrapper DoD
- [ ] Consumes T6 ME-branch-length tree (no per-job FastTree re-run); confirmed thread flag from T2

**Acceptance:** Wrapper DoD.

**Dependencies:** T8

**Agent notes:** GPL-3 + FastTree GPL-2+. APPLES (v1) vs APPLES-2: use APPLES-2 (`apples` 2.0.11).

---

### T15: DIPPER wrapper + Newick→placement adapter
**Milestone:** M1
**Description:** Run DIPPER `--add` in aligned and raw modes, CPU and GPU builds, and convert its Newick output into per-query placements.

**Deliverables:**
- [ ] Wrapper DoD for arms `dipper-m` (`-i m`, shared alignment) and `dipper-r` (`-i r`, raw)
- [ ] CPU build (bioconda `dipper`) and GPU build (Apptainer, pinned digest, `--gpu` measurement)
- [ ] `src/phylobench/dipper_adapter.py`: prune all other queries → attachment backbone edge → placement row
- [ ] Test: output tree restricted to backbone leaves has RF = 0 vs. input tree (fail loudly otherwise)
- [ ] Thread cap via cgroup/taskset (no native flag)

**Acceptance:** Wrapper DoD; adapter unit tests on hand-built trees; RF check passes on fixture.

**Dependencies:** T8

**Agent notes:** design §1.3. Document joint-placement caveat; do not "fix" it. MIT — the only fully non-copyleft candidate.

---

### M2 — Datasets (weeks 2–5, parallel with M1)

### T16: Similarity metric and query selection
**Milestone:** M2
**Description:** Compute per-query similarity to the backbone, bin it, and build query sets spanning the similarity gradient at fixed backbone size.

**Deliverables:**
- [ ] `src/phylobench/similarity.py`: `sim(q) = max_b identity(q,b)` over q-non-gap columns
- [ ] Bins: 75–80, 80–85, 85–90, 90–95, 95–97, 97–99, 99–99.9, ≥99.9
- [ ] Natural sampling + clade-pruning with backbone resubsampling to fixed N
- [ ] Per-bin lineage composition report

**Acceptance:** unit tests for identity/binning; on fixture, resubsampled backbone has exactly N leaves and no leaf ≥ hi identity to the query.

**Dependencies:** T5

**Agent notes:** design §4. Must be fast on 200k × 1k (vectorize; no O(n²) Python loops).

---

### T17: RNAsim dataset
**Milestone:** M2
**Description:** Download RNAsim and build backbones and query sets for primary accuracy and scaling.

**Deliverables:**
- [ ] Scripted download with checksum (`workflow/rules/data_rnasim.smk`)
- [ ] Subsamples 1k/10k/100k/200k; FastTree backbones on true alignment
- [ ] Query sets via T16 × 3 replicates; scaling query sets 1k/10k/100k on 10k backbone

**Acceptance:** dataset manifest (`data/rnasim/MANIFEST.tsv`) with counts per size/bin/replicate; bins with <100 queries flagged.

**Dependencies:** T16

**Agent notes:** Dryad doi:10.5061/dryad.78nf7dq — verify contents. Data in `data/` only, never committed.

---

### T18: GG2 full-length 16S dataset
**Milestone:** M2
**Description:** Build realistic 16S backbones from Greengenes2 with leave-one-out queries, full-length and in-silico V4 fragments.

**Deliverables:**
- [ ] Scripted GG2 backbone download (pinned release) + checksum
- [ ] Pruned GG2 topology subsamples (1k/10k/100k/max)
- [ ] LOO queries via T16; V4 fragments (length per Q4, default 150 nt) cut with 515F/806R primer coordinates

**Acceptance:** manifest per size/bin/query-type; LOO query never present in its own backbone (test).

**Dependencies:** T16

**Agent notes:** design §3.1. Use the genome-informed backbone, not DEPP-placed ASVs, for truth.

---

### T19: EMP V4 ASV query sets
**Milestone:** M2
**Description:** Assemble real EMP/Qiita V4 ASVs for realistic runtime/memory scaling and SEPP-concordance runs (no ground truth).

**Deliverables:**
- [ ] Scripted download of EMP Deblur ASVs (pinned source) + checksum
- [ ] Query sets 1k/10k/100k; per-ASV `sim` vs. GG2 backbone
- [ ] Evaluate rule outputs **concordance with SEPP** (edge distance to SEPP placement), labeled as such

**Acceptance:** manifest; concordance column never labeled "accuracy".

**Dependencies:** T18

**Agent notes:** GG2 ASV positions came from DEPP — not truth.

---

### T20: SARS-CoV-2 dataset
**Milestone:** M2
**Description:** Build a high-similarity stress-test panel from the UShER public MAT.

**Deliverables:**
- [ ] Pinned MAT date + checksum; `matUtils` subsample to 1k/10k/100k
- [ ] Reconstructed sequences from MAT + Wuhan-Hu-1; LOO queries
- [ ] Collapsed-polytomy truth sets for T5's edge-distance

**Acceptance:** manifest; ≥99% of queries fall in the ≥99% bins (sanity).

**Dependencies:** T16

**Agent notes:** design §3.1, Q8. Separate panel — never pooled with 16S.

---

### M3 — Runs (weeks 5–8)

### T21: HPC profile and pilot
**Milestone:** M3
**Description:** Create the SLURM execution profile and run a pilot at 1k/10k to calibrate resources and decide which tools attempt 100k/200k.

**Deliverables:**
- [ ] `profiles/slurm/config.yaml` (Snakemake SLURM executor): exclusive nodes, 16 cpus/task, per-tool mem/time, 256 GB / 48 h caps
- [ ] Verified cgroup v2 `memory.peak` on compute nodes (or documented fallback)
- [ ] Pilot results + `docs/pilot.md` with per-tool go/no-go for 100k/200k
- [ ] Fill `deployer.md` / `diagnostics.md` cluster specifics

**Acceptance:** pilot `runs.tsv` covers all arms at 1k/10k on RNAsim; timing repeat CV reported.

**Dependencies:** T9, T10, T11, T12, T13, T14, T15, T17

**Agent notes:** Q3 answers drive this. Never submit jobs without the user confirming the allocation/account.

---

### T22: Full benchmark execution
**Milestone:** M3
**Description:** Run all datasets × sizes × tools × arms on HPC and archive results.

**Deliverables:**
- [ ] `config/full.yaml`
- [ ] Complete `results/runs.tsv` + `results/accuracy.tsv`, DNFs recorded
- [ ] Archived results tarball with checksums and env locks

**Acceptance:** every planned (dataset, size, tool, arm, replicate) cell has a row with status; no silent gaps.

**Dependencies:** T21, T17, T18, T19, T20

**Agent notes:** Deployer owns submission; two-failed-attempt rule applies per job class.

---

### M4 — Analysis & report (weeks 8–10)

### T23: Benchmark plots
**Milestone:** M4
**Description:** Produce the final figures from the tidy result tables.

**Deliverables:**
- [ ] Accuracy (edge distance, delta error) vs. similarity bin, faceted by dataset × backbone size × query type
- [ ] Runtime and peak memory vs. similarity bin (per request), vs. backbone size, and vs. query count (log–log), DNF markers
- [ ] Placement-only (arm A) vs. end-to-end (arm B) comparison
- [ ] GPU memory panel (DIPPER)
- [ ] `workflow/scripts/plot_*.py`, reproducible from `results/*.tsv` alone

**Acceptance:** `snakemake plots` regenerates all figures from committed-sized sample TSVs in CI.

**Dependencies:** T22

**Agent notes:** Never pool bins across backbone sizes (design §4).

---

### T24: Report and recommendation
**Milestone:** M4
**Description:** Write the results report with SEPP non-inferiority per bin, cost table, license table, and a recommendation.

**Deliverables:**
- [ ] `results.md`: methods summary, per-bin non-inferiority vs. SEPP (Q9 rule), runtime/memory table, license table (from `LICENSES.md`), caveats (circular truth, DIPPER joint placement, confounds)
- [ ] Final recommendation with GPL/AGPL drawbacks called out

**Acceptance:** every claim cites a figure or table; user sign-off.

**Dependencies:** T23, T2

**Agent notes:** design §0 and §1.4.

---

## Done

<!-- Move tickets here once merged. Keep the T{N} prefix. -->
