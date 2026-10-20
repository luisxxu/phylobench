# phylobench — Design

**Status:** v1 — Planner ⇄ SWE **reconciled** (2 rounds); all user decisions resolved (§9), 2026-10-14.
**Goal:** Find a phylogenetic placement tool that is **at least as accurate as the lab's production SEPP** (production SEPP settings + GG2, run standalone), cheaper in **runtime and/or peak memory**, for eventual wide use across **Knight Lab infrastructure (Qiita / QIIME 2)**. Licenses are **reported, not used as a decision criterion**. Report accuracy, runtime, memory per tool across a **gradient of query-to-backbone similarity**.

Markers: **[V]** verified against upstream repo / GitHub API / bioconda / PyPI / GG2 or EMP FTP (2026-10-13/14). **[UNVERIFIED]** must be confirmed in T2 or T27. **⚑** deliberate deviation or open reconciliation point.

---

## 0. TL;DR

1. **Licenses are flagged, not decided on.** DIPPER is the only end-to-end MIT option; SCAMPP family = MIT wrapper over GPL-3 (pplacer) or AGPL-3 (EPA-ng) engines; APPLES-2, pplacer, SEPP are GPL-3 (§1.4). The recommendation is made on accuracy, runtime, and memory only.
2. **Baseline = production SEPP exactly:** `run-sepp.sh -A 1000 -P 5000` on the **GG2 2022.10 SEPP reference** (331,269 leaves) [V]. All other tools get the same backbone tree, alignment, and model numbers.
3. **Two backbones:** **Bfull\*** (production reference minus held-out queries) and **B10k\*** (10k subsample). Query-count scaling comes from real EMP ASVs.
4. **Truth = held-out clades + singletons in the GG2 tree** — one pruned backbone, one ref-prep. EMP ASVs have no truth (GG2 placed them with DEPP [V]) → labeled *concordance with SEPP*.
5. **Fair alignment:** challengers get a single-HMM `hmmalign` alignment (BSD, deployable); oracle alignment is a B10k\* diagnostic; SEPP aligns itself.
6. **DIPPER emits Newick with joint placement** → scored via an O(N) adapter (T25). jplace arms are converted to inserted trees with `guppy tog` (T25), as production does.
7. **Memory:** cgroup v2 `memory.peak` if available, else psutil process-tree RSS sum; source recorded per row.

---

## 1. Tool matrix

### 1.1 Arms and pins

| Arm | Upstream | Install (linux-64) | Pin | Own license | Engine executed |
|---|---|---|---|---|---|
| `sepp-prod` | github.com/smirarab/sepp `run-sepp.sh` with production flags (as q2-fragment-insertion calls it) | bioconda sepp/pplacer/hmmer at QIIME 2 2026.7 versions, own locked env (not QIIME 2) | 4.5.6; pplacer 1.1.alpha19 | GPL-3 | SEPP → pplacer alpha19 + HMMER; `guppy tog` inside run |
| `sepp-pinned` (conditional ⚑) | github.com/smirarab/sepp | bioconda `sepp` (uses conda pplacer/hmmer, not bundled [V]) | 4.5.6 [V] | GPL-3 | SEPP → shared pinned pplacer |
| `pplacer` | github.com/matsen/pplacer | bioconda | **1.1.alpha23** (alpha20 has OCaml-5 `-j` deadlock [V]; alpha22→23 is build-only) | GPL-3 [V] | itself |
| `pplacer-scampp` | github.com/ewedell/BSCAMPP `run_scampp.py --placement-method pplacer` [V] | PyPI `bscampp` | 1.0.8 [V] | MIT [V] | pplacer (GPL-3) |
| `epa-ng-scampp` | same, `run_scampp.py` | PyPI `bscampp` | 1.0.8 | MIT | EPA-ng (AGPL-3) |
| `bscampp` | same, `run_bscampp.py` (subtree 2000, votes 5 [V]) | PyPI `bscampp` | 1.0.8 | MIT | EPA-ng (AGPL-3) |
| `epa-ng` (control, B10k\* only) | github.com/pierrebarbera/epa-ng | bioconda | 0.3.8 [V] | AGPL-3 [V] | itself |
| `apples2` | github.com/balabanmetin/apples | PyPI `apples` | 2.0.11 [V] | GPL-3 [V] | itself (+ FastTree GPL-2+ for ME lengths) |
| `dipper-m` (CPU + GPU) | github.com/TurakhiaLab/DIPPER `--add -i m` | CPU: bioconda `dipper` (`dipper_cpu` only, issue #13 [V]); GPU: `swalia14/dipper` → Apptainer, digest-pinned | 0.1.5 [V] | MIT [V] | itself |
| `dipper-r` (B10k\* only) | same, `-i r` (MinHash k=15, sketch 1000 [V]) | same | 0.1.5 | MIT | itself |

- QIIME 2 2026.7 ships pplacer alpha19 ≠ our alpha23 → `sepp-pinned` required (T2, `docs/audit.md`).
- `sepp-pinned` runs only if T2 finds the QIIME 2 env's sepp/pplacer/hmmer versions differ from our pins (otherwise it duplicates `sepp-prod`). T2 records those versions and the binary paths SEPP resolves from `main.config` in `docs/audit.md` with the explicit rule "equal to pins → no sepp-pinned arm".
- All pplacer-based arms except `sepp-prod` share one pplacer binary.
- BSCAMPP/SCAMPP use external pinned `pplacer`/`epa-ng`/`taxit` and an immutable per-env config written at env deploy — never `~/.bscampp` or bundled `bscampp/tools/*` binaries [V]. (Exception: the bundled `tools/hamming_distance/*` helpers execute.)
- `ewedell/SCAMPP` (old, no license file [V]) is not used.

### 1.2 Inputs (all arms share Bfull\* or B10k\*)

| Arm | Query input | Parameter file |
|---|---|---|
| sepp-prod / sepp-pinned | unaligned | shipped `raxml-info.txt` (reused unchanged on pruned trees) |
| pplacer, pplacer-scampp | aligned (hmmalign; oracle at B10k\*) | shipped `raxml-info.txt` (what production pplacer reads via SEPP) |
| epa-ng, epa-ng-scampp, bscampp | aligned | `GTR{rates}+FU{freqs}+G4m{1.0}` converted from the info file (T6) ⚑ |
| apples2 | aligned | tree with ME branch lengths (FastTree `-nome -noml`, once per backbone) — required by APPLES [V] |
| dipper-m / dipper-r | aligned / raw | none (JC; CPU `dipper_cpu` defaults to `-d 1` despite its help text, so we pass `-d 2` explicitly) |

⚑ **Shipped model is GTRCAT with alpha=1.0 [V] — no fitted Γ.** Every likelihood arm gets the numbers production SEPP's pplacer gets. Deliberate deviation from per-family re-estimation; stated in the report. Alpha=1.0 handicaps *absolute* accuracy for all likelihood arms equally; the paired verdict vs. sepp-prod stays valid. **T6 must test** that the EPA-ng model string parses in raxml-ng/EPA-ng at the pinned versions and that the info-file rate order maps correctly to AC/AG/AT/CG/CT/GT — a silent permutation would corrupt every EPA-ng-family number.

### 1.3 DIPPER adapter and jplace→tree conversion (T25)

- **Adapter (Newick → placement), one O(N+Q) postorder pass:** assign each backbone leaf a random 64-bit hash; XOR-signature every backbone edge; in DIPPER's output (any rooting: edges are matched by canonical splits, so the adapter is rooting-agnostic) compute backbone-only signatures per node; each query walks up to the first ancestor with non-empty signature → backbone edge. Fail-loud RF check that the backbone topology is unchanged [UNVERIFIED whether `--add` preserves it]. Joint-placement caveat reported; delta error computed at B10k\*.
- **jplace → tree:** `guppy tog` from the pinned pplacer env, measured as phase `convert` for every jplace arm (`sepp-prod` runs it internally, as production does).

### 1.4 Licenses (flagged only — no license analysis, no decisions)

Per user decision, the report lists each arm's own license and the license of the engine it executes (table in §1.1, `LICENSES.md`) and flags copyleft (GPL/AGPL) as a note for later review. Licenses do not affect scope, ranking, or the recommendation.

---

## 2. Fairness and protocol

| Decision | Choice |
|---|---|
| Backbone | **Bfull\*** = shipped GG2 2022.10 SEPP reference minus held-out queries (~3% smaller); **B10k\*** = uniform 10k subsample of Bfull\*. Tree and alignment pruned, never re-estimated; pruning suppresses degree-2 nodes and keeps branch lengths additive. The pruned reference's checksum is recorded in the dataset manifest so all arms provably share it. ⚑ |
| Model params | Shipped numbers for every likelihood arm (§1.2); ME branch lengths for APPLES-2 only. Ref-prep timed as a one-time cost. |
| Alignment | **Primary:** `hmmalign` against one `hmmbuild` profile of the backbone alignment, masked to backbone columns (the realistic deployable path); shared by every aligned-input arm (phase `align`, timed). **Diagnostic (B10k\* only):** oracle (region's own GG2 alignment columns). SEPP aligns internally. A poor single-profile alignment lowers challengers vs. SEPP — reported as a property of the alignment arm, not the tools; the oracle-vs-hmmalign gap at B10k\* is the evidence. ⚑ |
| End-to-end cost | `sepp-prod` total vs. `align + place + convert` for each challenger. |
| Regions | full-length, **V4-150**, **V3–V4**, **V1–V2**, **V1–V3**, **V4–V5**, **V5–V6**, **V6–V8**, **V7–V9**, cut by primer matching on ungapped GG2 sequences (V4-150: 515F Parada `GTGYCAGCMGCCGCGGTAA` / 806R Apprill `GGACTACNVGGGTWTCTAAT`, first 150 nt after 515F; V3–V4: 341F `CCTACGGGNGGCWGCAG` / 805R `GACTACHVGGGTATCTAATCC`; V1–V2: 27F `AGAGTTTGATCMTGGCTCAG` / 338R `GCTGCCTCCCGTAGGAGT`; V1–V3: 27F / 534R `ATTACCGCGGCTGCTGG`; V4–V5: 515F Parada / 926R Parada `CCGYCAATTYMTTTRAGTTT`; V5–V6: 784F `AGGATTAGATACCCTGGTA` / 1064R `CGACRRCCATGCANCACCT`; V6–V8: 926F `AAACTYAAAKGAATTGRCGG` / 1392R `ACGGGCGGTGTGTRC`; V7–V9: 1115F `CAACGAGCGCAACCCT` / 1492R `GGTTACCTTGTTACGACTT`; exact match, primers excluded from the amplicon) and projected onto alignment columns. Sequences without primer hits excluded and counted. |
| Threads | 16/job (SLURM + each tool's flag). DIPPER native `-T {threads}` at 0.1.5 (TBB global_control) is CPU only; the GPU build has no thread flag, so the harness records 16 and SLURM bounds it. BSCAMPP `--threads 16 --cpus-per-job 2`. |
| SEPP | `run-sepp.sh <q> <name> -x 16 -A 1000 -P 5000 -a -t -r` [V], one working dir and `-p`/`-d` temp dirs per job. |
| Placement scored | Top-LWR only. |

---

## 3. Truth and metrics

| Data | Role | Truth label |
|---|---|---|
| GG2 held-out clades + singletons (§4), 9 regions, 1,000 source sequences each (paired across regions) | **Accuracy** at B10k\* and Bfull\* | Attachment edge in GG2 tree restricted to backbone ∪ {q}; held-out-clade members share its stem edge. **Caveat:** GG2 16S were inserted by uDance (APPLES-2 + ML refinement) [V] → possible bias toward apples2 (and ML arms). An apples2 win is **not** read as a placement-quality win without the non-uDance cross-check (T3 simulator scaled to a 10k true-tree panel, in T16) or T17-lite. |
| EMP Deblur V4 150 bp ASVs (release1 [V]), 10k and 100k | **Runtime, memory, query-count scaling, concordance** | *Concordance with sepp-prod* (edge distance). Never "accuracy". |
| (optional) T17-lite RNAsim 10k | Bias-free truth | True tree (§9.8) |

**Metrics:** edge distance (primary; polytomy-aware: collapse zero-length branches, any edge in the true equivalence set scores 0); fraction exactly correct; delta error (B10k\* only); DNF recorded as data.

Edge distance is counted in the backbone with zero-length edges contracted: a true position that is a vertex (query attaches at a multifurcation) is at distance 0 from every incident edge, and two edges sharing a contracted vertex are adjacent (distance 1). Unplaced queries are NaN rows; consumers never `skipna` silently (`phylobench.accuracy`).

**Verdict ("as accurate as SEPP"):** per (region, backbone size), pooled across bins, paired per-query difference in edge distance (tool − sepp-prod); **pass** if the upper 95% bootstrap CI of the mean ≤ **0.25 edges**. Per-bin results are descriptive only (avoids the multiple-testing trap).

**Short-window caveat:** identical V4-150 fragments can fit many edges equally; edge distance then reflects unavoidable ambiguity — hence the *paired* comparison against SEPP.

---

## 4. Percent similarity and held-out design

**Definition:** `sim(q)` = max identity to any backbone leaf over q's non-gap columns in the primary alignment (within the region window for regions), computed **after** the held-out backbone is built, per backbone size. **Never pool bins across sizes.**

**Bins:** 75–80, 80–85, 85–90, 90–95, 95–97, 97–99, 99–99.9, ≥99.9.

**Held-out construction (one pass, one pruned backbone per size):**
1. Using GG2 backbone taxonomy (`2022.10.backbone.tax.qza` [V]), remove disjoint clades from the shipped tree + alignment: ~6 genera, 6 families, 4 orders, 3 classes, 2 small phyla, each ≤2k leaves. Higher ranks → lower-similarity queries.
2. Remove **~500 singleton leaves** outside those clades. In a 331k backbone singletons usually keep a near-identical neighbor (operon copies, strains, sister species) — this fills **99–99.9 and ≥99.9 on 16S** without SARS-CoV-2. If ≥99.9 at Bfull\* has <100 queries, replace trailing singletons with leaves whose nearest remaining neighbor is ≥99.9%, keeping the query count at 1,000 (`docs/gg2.md`). V4-150 fragments push identity further up.
3. Queries = 1,000 held-out sequences (~500 singletons + ~500 clade members, ≤25 per clade), same sources for every region.
4. **Bfull\*** ≈ 320k leaves (exact duplicate aligned rows removed after hold-out — FastTree needs none); **B10k\*** = uniform subsample of Bfull\*. At B10k\* the same queries shift to lower bins (sparser neighborhoods) — expected, reported per size.

**Confounds reported:** per-bin lineage composition; per-bin query length (regions); bins with <100 queries flagged, not padded.

---

## 5. Scale

| Axis | Values |
|---|---|
| Backbones | B10k\*, Bfull\* (≈320k) |
| Accuracy queries | 1,000 per region × 9 regions |
| Scaling queries | EMP V4 10k and 100k on both backbones |
| Arms at Bfull\* | sepp-prod, (sepp-pinned), pplacer-scampp, epa-ng-scampp, bscampp, apples2, dipper-m CPU + GPU, plus **one capped pplacer attempt** (1k V4 queries, 24 h cap enforced by the harness as `status=timeout`, fixed before the run, no retry) ⚑ |
| B10k\*-only arms | pplacer, epa-ng, dipper-r, oracle-alignment diagnostics |
| Replicates | One query set; bootstrap over queries. Timing CV from 3 repeats of one B10k\* config. |
| Limits | Exclusive x86 node, 16 threads, >256 GB node memory, 72 h wall. GPU: RTX 3090 24 GB. |

**Expected walls** [UNVERIFIED until pilot]: pplacer and plain EPA-ng at Bfull\*; DIPPER GPU at Bfull\* may exceed 24 GB.

---

## 6. Measurement

| Quantity | How |
|---|---|
| Wall / CPU time | Harness wall clock; CPU user+sys from `os.wait4` rusage |
| Peak memory | cgroup v2 `memory.peak` if T27 confirms; else psutil process-tree RSS **sum** sampled every 0.5 s; `os.wait4` max RSS (`ru_maxrss`, same as GNU `time -v`) always recorded; both recorded; `mem_source` column |
| GPU memory | `nvidia-smi --query-compute-apps=used_memory` every 0.5 s, max |
| Phases | `ref-prep`, `align`, `place`, `convert`, `evaluate` as separate rules, each measured; `sepp-prod` = one `align+place+convert` phase |
| Provenance | node, CPU model, GPU model, threads, version, env-lock hash, command |

Tests skip explicitly when cgroups are absent — never pass vacuously.

- `results/runs.tsv`: `dataset, backbone, region, query_n, tool, arm, aln, phase, wall_s, cpu_s, peak_mem_mb, mem_source, gpu_peak_mb, status, threads, node, cpu_model, version, run_id, dnf_reason, env_lock_sha256`; `status` ∈ {`ok`, `DNF`} (`dnf_reason` = the harness status: `oom|timeout|error`)
- `results/accuracy.tsv`: `run_id, query, sim, bin, edge_dist, exact, delta_err, label(accuracy|concordance)`

---

## 7. Stack

| Concern | Choice |
|---|---|
| Workflow | Snakemake; Snakefile includes `rules/*.smk` and loads `config/arms/*.yaml` by glob (pre-registered in T8) so wrappers only add files |
| Environments | One pinned env + lockfile per tool in `envs/`; `sepp-prod` pins QIIME 2 2026.7's sepp/pplacer/hmmer/openjdk versions in its own env, without QIIME 2; DIPPER GPU as digest-pinned Apptainer `.sif` |
| Code | `src/phylobench/` pure modules: `measure`, `jplace`, `accuracy`, `similarity`, `regions`, `heldout`, `dipper_adapter`, `models`. Snakemake rules do I/O. |
| Gates | Local (darwin agents): `pytest`, `ruff check`, `ruff format --check`, `snakemake --lint`. CI (`ubuntu-latest`): same + smoke test with every arm (`dipper_cpu`). Unit tests need no tool binaries. |
| Compute | SLURM x86 exclusive nodes; one RTX 3090 GPU node type; job submission needs explicit user confirmation |
| Data | `data/`, `results/`, `.snakemake/` gitignored; scripted checksummed downloads (GG2 2022.10, EMP release1) |

---

## 8. Timeline — 5 weeks (3-week floor below)

| Week | Work |
|---|---|
| W1 | **T27** day-0 cluster check (human). T1, T2, T3, T4, T5, T16, T25 |
| W2 | T6, T7, T8 (smoke). T18 (Bfull\*/B10k\*, queries). **T26 sepp-prod** first; T10, T12, T14, T15 start |
| W3 | Wrappers merge. T19 (EMP). **T21 pilot** (DIPPER GPU, sepp-prod first, others as they land). T23 plot scripts on sample TSVs |
| W4 | T22 staged runs: B10k\* accuracy → Bfull\* accuracy → EMP scaling. Optional T17-lite only if B10k\* is done by W4 day 2 |
| W5 | Buffer; T23 final figures; T24 report |

**Critical path:** T1 → T2 → T6/T7 → T8 → T26 → T21 → T22 → T23 → T24 (T27 before T21).

**3-week floor (approved) — cut in order:** T17-lite → `dipper-r` → oracle diagnostic → the six non-V4 amplicon regions (keep full-length + V4-150) → EMP 100k → `sepp-pinned`. Runs start W3; report end of W3.

---

## 9. Resolved decisions (user, 2026-10-13/14)

1. **Licenses** are flagged in the report only; no license analysis and no decisions based on license (§1.4).
2. **Production SEPP** = q2-fragment-insertion + GG2 → `sepp-prod` on the 2022.10 SEPP reference (only GG2 SEPP reference; 2024.09 has none [V]); `-A 1000 -P 5000` [V]; env = **latest QIIME 2 amplicon distribution**, sepp/pplacer/hmmer versions recorded (T2). Amended 2026-10-19 (user): standalone run-sepp.sh, not q2-fragment-insertion/QIIME 2 env; versions mirror QIIME 2 2026.7.
3. **Compute:** SLURM x86, ≥16 cores, >256 GB/node, 72 h; RTX 3090 24 GB. cgroup support unknown → T27 + psutil fallback.
4. **Regions:** full-length, V4-150, V3–V4, V1–V2, V1–V3, V4–V5, V5–V6, V6–V8, V7–V9 (the last five added 2026-10-16).
5. **Defaults accepted:** 10k-ASV scaling (to 100k); GG2 backbone; EPA-ng control; 0.25-edge margin (per region × size, pooled); DEPP out of scope; GPU shown separately, CPU ranks; inserted tree needed downstream (→ T25 via `guppy tog`).
6. **Deadline:** 3–5 weeks → 5-week plan; the **3-week floor cut order (§8) is approved**.
7. **SARS-CoV-2:** deferred out of this cycle (T20); high-similarity bins are filled on 16S (§4).
8. **RNAsim:** optional **T17-lite** (one 10k panel, accuracy only) in W4 if on schedule; first cut otherwise. APPLES-2 is not expected to be a strong contender; the T16 simulated cross-check still runs.

---

## 10. Risks

| Risk | Mitigation |
|---|---|
| 5-week schedule; cluster surprises (cgroups, Apptainer, CUDA on 3090) | T27 day 0; staged runs; cut list |
| Pruned Bfull\* ≠ shipped production reference | ~3% smaller; stated |
| Shipped model GTRCAT alpha=1.0 (no fitted Γ) | Same numbers for all arms (matches production); stated |
| GG2 truth seeded by APPLES-2 (uDance) | Report caveat; 10k simulated true-tree cross-check (T16); optional RNAsim-lite |
| Single-HMM alignment poor on V1–V2 / divergent clades | Penalizes all challengers equally; oracle diagnostic separates alignment- from placement-limited error |
| DIPPER `--add` alters backbone; GPU memory at 320k; v0.1.x churn | RF check; digest pin; CPU fallback |
| BSCAMPP bundled binaries / `~/.bscampp` | Immutable per-env config written at env deploy; provenance shows binary paths |
| psutil fallback misses short spikes | `mem_source` column; caveat on memory panels |
| Primer matching drops divergent-primer sequences | Count and report exclusions per region |
| EMP biom observation IDs may not be sequences [UNVERIFIED] | T19 verifies first |
