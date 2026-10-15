# Tickets

One ticket per PR. Move completed tickets to **## Done**. Ticket `T{N}` ↔ GitHub issue `#{N}` in the `phylobench` Project. Design rationale lives in [`design.md`](design.md) (v1) — read the cited sections before planning.

Numbers are never reused or renumbered. Merged or deferred tickets keep their number under **## Merged / Deferred** (issue labeled `deferred` or closed as merged).

Each ticket has: **Week**, **Description**, **Deliverables**, **Acceptance** (what QA checks), **Dependencies**, **Agent notes**.

**Wrapper DoD** (T10, T12, T14, T15, T26):
- [ ] `envs/<tool>.yaml` pinned per `design.md` §1.1 + lockfile (linux-64)
- [ ] `workflow/rules/<tool>.smk` + `config/arms/<arm>.yaml` — **add files only**; never edit shared config (T8 globs them)
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

### T2: License and install audit
**Week:** W1
**Description:** Install every in-scope tool from a pinned env on linux-64 (in CI), record licenses (flag only), and resolve the install-level `[UNVERIFIED]` items.

**Deliverables:**
- [ ] `envs/{pplacer,sepp,bscampp,apples,dipper,epa-ng,fasttree,hmmer}.yaml` + lockfiles; `envs/qiime2.yaml` locked to the latest QIIME 2 amplicon distribution
- [ ] CI job printing each tool's version
- [ ] `LICENSES.md` complete up front (wrappers never edit it): tool → own license → engine executed → copyleft flag → source URL at pinned tag. No license analysis.
- [ ] `docs/audit.md`: QIIME 2 env's sepp/pplacer/hmmer versions and the binary paths SEPP resolves from `main.config`, with the rule **"equal to our pins → no `sepp-pinned` arm"**; pplacer pin (alpha22) confirmation; proof BSCAMPP/SCAMPP use the external pinned binaries + per-job config; APPLES-2 thread flag; 

**Acceptance:** every env solves and prints its pinned version in CI; `LICENSES.md` has a row for every arm in design §1.1 plus HMMER, FastTree, guppy.

**Dependencies:** T1

**Agent notes:** Read LICENSE files at the pinned tag; record what *executes*. Licenses are flagged only — never a reason to change scope. No cluster/GPU work here (that's T27).

---

### T3: Tiny fixture dataset
**Week:** W1
**Description:** Pure-Python seeded simulator producing a tiny dataset with known truth, regions, and held-out clades for unit and smoke tests.

**Deliverables:**
- [ ] `workflow/scripts/make_fixture.py` (numpy GTR simulation on a random ~100-leaf tree, gap-free, ~1,500 sites, fixed seed, fixed-precision output)
- [ ] `tests/fixtures/tiny/`: true tree, alignment, primer sites embedded so all four regions extract, held-out clades + singletons, `truth.tsv`, a GTRCAT-style `raxml-info.txt`
- [ ] Simulator importable so T16 can scale it to a 10k-leaf cross-check panel

**Acceptance:** regeneration reproduces committed files byte-identically; total < 1 MB.

**Dependencies:** T1

**Agent notes:** No downloads. Mirror the GG2 SEPP reference file shapes so downstream code is exercised realistically.

---

### T4: Measurement harness
**Week:** W1
**Description:** One wrapper that runs a command and records wall, CPU, peak memory (with source), GPU memory, status, and provenance as JSON.

**Deliverables:**
- [ ] `src/phylobench/measure.py` + CLI `phylobench-measure --out run.json [--gpu] [--timeout S] -- <cmd>`
- [ ] Peak memory: cgroup v2 `memory.peak` if readable for the job's own cgroup, else psutil process-tree RSS **sum** sampled every 0.5 s; `/usr/bin/time -v` max RSS always recorded; `mem_source` field
- [ ] `--timeout` enforced by the harness → `status=timeout` (used for the capped pplacer Bfull\* attempt)
- [ ] Status: `ok | oom | timeout | error`; provenance: node, CPU/GPU model, threads, version, env-lock hash, command

**Acceptance:** tests for status classification and JSON schema; a forking command reports peak ≥ sum of children via psutil; cgroup tests `pytest.skip` explicitly when unavailable (never pass vacuously).

**Dependencies:** T1

**Agent notes:** design §6. `time -v` reports the largest single process — never the silent primary.

---

### T5: Accuracy metrics library
**Week:** W1
**Description:** Pure functions to read placements and score them against a true tree.

**Deliverables:**
- [ ] `src/phylobench/jplace.py`: minimal reader → top-LWR (query, edge_num, distal, pendant)
- [ ] `src/phylobench/accuracy.py`: polytomy-aware edge distance (collapse zero-length branches; any edge in the true equivalence set = 0), exact-correct flag, delta error (FN bipartitions)
- [ ] Unit tests on hand-built trees incl. polytomies, zero-length edges, degree-2 suppression

**Acceptance:** tests pass; edge distance for 1k queries on a 330k-leaf random tree runs in < 60 s.

**Dependencies:** T1

**Agent notes:** design §3. No O(n²) loops over leaves.

---

### T16: Similarity, held-out construction, and simulated cross-check
**Week:** W1–W2
**Description:** Build the held-out backbone + query sets, compute per-query similarity bins, and produce the non-uDance simulated cross-check panel.

**Deliverables:**
- [ ] `src/phylobench/heldout.py`: choose disjoint clades at mixed ranks (~6 genera, 6 families, 4 orders, 3 classes, 2 small phyla, each ≤2k leaves) + ~500 singletons; ≥99.9 seeding rule; prune tree (suppress degree-2 nodes, additive branch lengths) and alignment once
- [ ] `src/phylobench/similarity.py`: `sim(q)` over q's non-gap columns (region window aware); bins 75–80 … ≥99.9; chunked numpy
- [ ] Per-bin lineage composition and query-length report
- [ ] 10k-leaf simulated true-tree panel from the T3 simulator (cross-check for uDance/APPLES-2 bias)

**Acceptance:** on the fixture, no held-out leaf remains in the backbone, every query's truth edge exists in the pruned tree, no degree-2 nodes remain; `sim` on 10k × 10k completes in < 5 min.

**Dependencies:** T3, T5

**Agent notes:** design §4. One pruned backbone per size — never per query.

---

### T25: Tree ↔ placement conversion (guppy tog + DIPPER adapter)
**Week:** W1–W2
**Description:** Convert jplace arms to inserted trees with `guppy tog` (as production does), and write the O(N) DIPPER Newick→placement adapter.

**Deliverables:**
- [ ] `workflow/rules/convert.smk`: `guppy tog` from the pinned pplacer env for every jplace arm (phase `convert`, harness-wrapped)
- [ ] `src/phylobench/dipper_adapter.py`: random 64-bit leaf hashes → XOR edge signatures → one postorder pass maps each query to its backbone edge; RF check that the backbone topology is unchanged (fail loudly)

**Acceptance:** `guppy tog` produces an inserted tree containing every query for a fixture jplace (CI); adapter maps 100k queries on a synthetic 330k tree in < 60 s; adapter unit tests on hand-built trees.

**Dependencies:** T2, T5

**Agent notes:** design §1.3. Document DIPPER joint placement; don't "fix" it. The adapter is pure Python and can start before T2; the convert rule needs T2's pplacer env.

---

### Week 2 — Ref-prep, alignment, smoke test, datasets

### T6: Reference prep (pruning + model conversion)
**Week:** W2
**Description:** Produce each tool family's parameter inputs from the shipped GG2 SEPP reference numbers — no re-estimation.

**Deliverables:**
- [ ] `src/phylobench/models.py`: parse `raxml-info.txt` (GTR rates, base freqs, alpha) → EPA-ng/raxml-ng string `GTR{rates}+FU{freqs}+G4m{1.0}`
- [ ] FastTree `-nome -noml` ME branch lengths for APPLES-2, once per backbone
- [ ] Rules wrapped by harness, phase `ref-prep`

**Acceptance:** model string parses in pinned raxml-ng and EPA-ng; test that info-file rate order maps to AC/AG/AT/CG/CT/GT (a known permuted input fails); pplacer and epa-ng each complete on the fixture using T6 outputs (CI).

**Dependencies:** T2, T3

**Agent notes:** design §1.2 ⚑ — alpha=1.0 is as shipped; don't "fix" it.

---

### T7: Regions and query alignment
**Week:** W2
**Description:** Extract amplicon regions by primer matching and produce the shared hmmalign alignment plus the oracle diagnostic.

**Deliverables:**
- [ ] `src/phylobench/regions.py`: primer matching (IUPAC-aware) on ungapped sequences for V4-150 (515F/806R, trimmed to 150), V3–V4 (341F/805R), V1–V2 (27F/338R) [confirm primer set], projected onto the row's alignment columns; exclusion counts
- [ ] Check aligned rows ungap to the full-length FASTA
- [ ] `hmmbuild` profile of each backbone + `hmmalign` + mask to backbone columns (phase `align`, timed)
- [ ] Oracle passthrough (region's GG2 alignment columns) for the B10k\* diagnostic

**Acceptance:** on the fixture, every region extracts for every sequence; aligned queries have exactly the backbone's column count.

**Dependencies:** T2, T3

**Agent notes:** design §2. WITCH is out of scope.

---

### T8: Tiny end-to-end smoke test
**Week:** W2
**Description:** Wire fixture → ref-prep → align → pplacer → convert → evaluate → tidy TSVs in CI, and set the add-files-only convention for wrappers.

**Deliverables:**
- [ ] `config/smoke.yaml` + `config/arms/pplacer.yaml`; minimal `workflow/rules/pplacer.smk`
- [ ] `workflow/rules/evaluate.smk` → `results/runs.tsv` and `results/accuracy.tsv` (columns per design §6)
- [ ] CI job `smoke`: `snakemake --use-conda --configfile config/smoke.yaml -c 2`, matrix driven by the `config/arms/` glob

**Acceptance:** CI smoke green; `accuracy.tsv` has one row per fixture query × region; pplacer full-length exact-correct ≥ 8/10 on the fixture.

**Dependencies:** T3, T4, T5, T6, T7

**Agent notes:** The template every wrapper copies. Adding a wrapper must not require editing any existing file.

---

### T18: GG2 16S dataset
**Week:** W2
**Description:** Download the GG2 2022.10 SEPP reference + taxonomy and build Bfull\*, B10k\*, and the held-out query sets for four regions.

**Deliverables:**
- [ ] Scripted checksummed download of `2022.10.backbone.sepp-reference.qza` and `2022.10.backbone.tax.qza`; unzip `.qza` without needing QIIME 2
- [ ] Bfull\* (shipped minus held-out, via T16) and B10k\* (uniform subsample of Bfull\*); **pruned-reference checksums recorded** in `data/gg2/MANIFEST.tsv`
- [ ] 1,000 held-out source sequences × {full-length, V4-150, V3–V4, V1–V2}; per-size bin counts; ≥99.9 seeding applied

**Acceptance:** manifest lists counts per backbone × region × bin and checksums; no held-out leaf in either backbone (test); bins < 100 flagged.

**Dependencies:** T7, T16

**Agent notes:** design §3–§4. `data/` never committed.

---

### T26: sepp-prod wrapper (production baseline)
**Week:** W2 (first wrapper)
**Description:** Run SEPP exactly as q2-fragment-insertion does in production — the yardstick for every verdict.

**Deliverables:**
- [ ] Wrapper DoD (env = locked QIIME 2 distribution — the one pinned-env exception)
- [ ] `run-sepp.sh <q> <name> -x {threads} -A 1000 -P 5000 -a -t -r` directly, per-job working dir and `-p`/`-d` temp dirs; single phase `align+place+convert`
- [ ] q2-vs-direct equality check: on a small query set against the **as-shipped** reference, `qiime fragment-insertion sepp` and the direct call give identical placements

**Acceptance:** Wrapper DoD; equality check passes in CI on the fixture-sized reference (and once on the shipped reference during T21).

**Dependencies:** T2, T8

**Agent notes:** design §2 SEPP row. Don't force our pinned pplacer here — production is what it is.

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

**Agent notes:** alpha22 pin (alpha20 deadlocks with `-j`).

---

### T12: bscampp family + EPA-ng control
**Week:** W2–W3
**Description:** `pplacer-scampp`, `epa-ng-scampp`, `bscampp`, and plain `epa-ng` from one pinned `bscampp` install.

**Deliverables:**
- [ ] Wrapper DoD for each of the four arms
- [ ] External pinned `pplacer`/`epa-ng` (never bundled `bscampp/tools/*`), per-job config dir (never `~/.bscampp`)
- [ ] `bscampp`: `--threads {threads} --cpus-per-job 2`, subtree 2000, votes 5 recorded explicitly
- [ ] `epa-ng` control configured for B10k\* only

**Acceptance:** Wrapper DoD; provenance shows external binary paths; two concurrent smoke runs don't share config.

**Dependencies:** T8

**Agent notes:** MIT wrapper, GPL-3 / AGPL-3 engines — `docs/tools/` must say which.

---

### T14: APPLES-2 wrapper
**Week:** W2–W3
**Description:** APPLES-2 on hmmalign-aligned queries against the ME-branch-length backbone.

**Deliverables:**
- [ ] Wrapper DoD; consumes T6 ME tree (no per-job FastTree)

**Acceptance:** Wrapper DoD.

**Dependencies:** T8

**Agent notes:** Results carry the uDance-bias caveat. Not expected to be a strong contender — keep the wrapper minimal.

---

### T15: DIPPER wrapper
**Week:** W2–W3
**Description:** DIPPER `--add` in aligned (`-i m`) and raw (`-i r`, B10k\* only) modes, CPU and GPU builds.

**Deliverables:**
- [ ] Wrapper DoD for `dipper-m` (CPU + GPU) and `dipper-r`
- [ ] CPU: bioconda `dipper`; GPU: Apptainer `.sif` pinned by digest, `--nv`, harness `--gpu`
- [ ] Output through T25 adapter; thread cap via cgroup/taskset

**Acceptance:** Wrapper DoD (CPU in CI); adapter RF check passes on fixture output.

**Dependencies:** T8, T25

**Agent notes:** First in the T21 pilot (GPU is new for the lab). GPU path can only be verified on the cluster (manual tier).

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
- [ ] Pilot at B10k\*: dipper-m (GPU, CPU) and sepp-prod first, then each merged arm; q2-vs-direct check on the shipped reference
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

**Acceptance:** every planned (backbone, region, arm) cell has a row with status — no silent gaps.

**Dependencies:** T21, T18, T19, T10, T12, T14, T15, T26

**Agent notes:** Deployer owns submission; two-failed-attempt rule per job class. Optional T17-lite only if B10k\* is done by W4 day 2.

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

### T1: Repo skeleton and CI — merged in #29 (2026-10-15)
Skeleton, `config/schema.yaml`, `envs/phylobench.yaml`, CI job `lint-test` (now a required check on `main`). Pin change: pandas 2.3.3 (bioconda snakemake 9.27.0 requires pandas <3).
