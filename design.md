# phylobench — Design

**Status:** Draft v0 (Planner). Open questions in §9 carry defaults so work can proceed; answers override defaults.
**Goal:** Find a phylogenetic placement tool that is **at least as accurate as SEPP**, improves on **runtime and/or peak memory**, and is ideally **not GPL**. Report accuracy, runtime, and memory per tool across a **gradient of query-to-backbone similarity**.

Verification markers: **[V]** = verified against upstream repo / GitHub API / bioconda / PyPI on 2026-10-13. **[UNVERIFIED]** = from memory or inference; must be confirmed in T2.

---

## 0. TL;DR — opinions to push back on

1. **Strictly, none of the 7 tools is non-GPL except DIPPER.** The SCAMPP family is MIT, but executes (and bundles) pplacer (GPL-3) or EPA-ng (**AGPL-3**). APPLES-2 is GPL-3. We need to know *why* GPL matters (Q1) before license drives the decision.
2. **DIPPER does not emit jplace.** It outputs a Newick tree with queries grafted in and places queries *jointly*. It needs its own adapter, and its accuracy must be scored carefully (§1.3).
3. **"Truth" is circular for two candidate datasets.** EMP V4 ASVs in Greengenes2 were *placed by DEPP* [V]; the SARS-CoV-2 public MAT was *built by UShER placement*. Primary accuracy should come from **RNAsim (true tree)** + **leave-one-out on full-length 16S backbones**; real EMP ASVs are for runtime/memory and SEPP-concordance.
4. **Runtime/memory barely depend on similarity** — they depend on backbone size and query count. Recommendation: plot **accuracy vs. similarity**, and **runtime/memory vs. backbone size and query count**. (A runtime-vs-similarity panel is still produced for completeness, but expect flat lines.)
5. **`/usr/bin/time -v` under-reports memory for tools that fork workers** (pplacer `-j`, SEPP, BSCAMPP). Use cgroup v2 `memory.peak` as primary.
6. **Add plain EPA-ng as an 8th "control" arm.** Nearly free, and the only way to separate "EPA-ng engine" from "SCAMPP/BSCAMPP framework" effects (Q7).

---

## 1. Tool matrix

### 1.1 Summary

| # | Tool | Upstream | Install (linux-64) | Pin | License | Engine it executes |
|---|---|---|---|---|---|---|
| 1 | pplacer | github.com/matsen/pplacer | bioconda `pplacer` | **1.1.alpha22** | GPL-3.0 [V] | itself |
| 2 | SEPP | github.com/smirarab/sepp | bioconda `sepp` | **4.5.6** [V] | GPL-3.0 [V] | pplacer + HMMER (BSD) |
| 3 | pplacer-SCAMPP | github.com/ewedell/BSCAMPP (`run_scampp.py --placement-method pplacer`) [V] | PyPI `bscampp` (not on bioconda [V]) | **1.0.8** [V] | MIT [V] | pplacer (GPL-3) |
| 4 | EPA-ng-SCAMPP | same repo, `run_scampp.py` (EPA-ng default) [V] | PyPI `bscampp` | 1.0.8 | MIT [V] | **EPA-ng (AGPL-3)** |
| 5 | BSCAMPP | same repo, `run_bscampp.py` [V] | PyPI `bscampp` | 1.0.8 | MIT [V] | EPA-ng (AGPL-3) default; pplacer optional |
| 6 | APPLES-2 | github.com/balabanmetin/apples | PyPI `apples` (not on bioconda [V]) | **2.0.11** [V] | GPL-3.0 [V] | itself; runs **FastTree (GPL-2+)** for ME branch lengths [V] |
| 7 | DIPPER | github.com/TurakhiaLab/DIPPER | bioconda `dipper` (CPU build only, issue #13 [V]); GPU via Docker `swalia14/dipper` → Apptainer | **0.1.5** [V] | MIT [V] | itself (Boost, TBB, CUDA runtime) |
| (8) | EPA-ng (control) | github.com/pierrebarbera/epa-ng | bioconda `epa-ng` | **0.3.8** [V] (not 0.3.9-pre) | AGPL-3.0 [V] | itself |

Notes:
- `ewedell/SCAMPP` is the old evaluation-scripts repo with **no license file** [V]. SCAMPP now lives inside `ewedell/BSCAMPP` (originally from `chry04/PLUSplacer`, MIT [V]). Do not use the old repo.
- **pplacer pin:** bioconda has alpha17/19/20/22/23 [V]. alpha20 moved to OCaml 5; alpha21 fixed a multiprocessing deadlock under OCaml 5 [V] → alpha20 unsafe with `-j`. alpha23 (2026-10-07) is too fresh. Pin **alpha22** and force SEPP + both pplacer-based SCAMPP variants to use **the same binary**. T2 decides whether alpha19 (used in published SCAMPP/SEPP results) is also needed for literature comparability.
- **BSCAMPP ships bundled `epa-ng`/`pplacer` binaries** under `bscampp/tools/` [V] and writes `~/.bscampp/main.config` on first run [V]. Risks: silently using bundled (unpinned) binaries; concurrent jobs racing on shared `$HOME`. Wrapper must set an explicit per-job config and binary paths.

### 1.2 Inputs, hardware, outputs

| Tool | Backbone aln | Backbone tree | Model/params file | Query input | Threads | GPU | Output |
|---|---|---|---|---|---|---|---|
| pplacer | yes | yes | FastTree log or RAxML-7/8 info / taxtastic refpkg; does not read raxml-ng bestModel [UNVERIFIED] | aligned, backbone width | `-j` (forks) | no | jplace |
| SEPP | yes | yes | RAxML info (`-r`) [V]; FastTree log? [UNVERIFIED] | **unaligned** (internal HMM ensemble) [V] | `-x` | no | jplace + ext. alignment |
| pplacer-SCAMPP | yes | yes | FastTree/RAxML log [V] | aligned (`-q`) [V] | `--threads` | no | jplace |
| EPA-ng-SCAMPP | yes | yes | raxml-ng `.bestModel` [V] | aligned | `--threads` | no | jplace |
| BSCAMPP | yes | yes | raxml-ng `.bestModel` [V] | aligned | `--threads`, `--cpus-per-job` (default 2) [V] | no | jplace |
| APPLES-2 | yes | yes, **ME branch lengths** (FastTree `-nome -noml`, default) [V] | none (JC69) | aligned or distance matrix [V] | `-T` [UNVERIFIED] | no | jplace [V] |
| DIPPER `--add` | one FASTA of backbone **+** queries, aligned (`-i m`) or raw (`-i r`, MinHash k=15, sketch 1000) [V] | yes (`-t`) [V] | none (`-d`, JC default) [V] | aligned or unaligned | **no documented thread flag** [V] → cap via cgroup/taskset | optional (CUDA/HIP); `dipper_cpu` [V] | **Newick only** [V] |
| EPA-ng | yes | yes | raxml-ng model/bestModel [V] | aligned | `-T` / `OMP_NUM_THREADS` [V] | no | jplace |

### 1.3 DIPPER challenges

- **No jplace.** Adapter (T15): for each query, prune all *other* queries, find the backbone edge the query hangs from, emit a jplace-equivalent row.
- **Joint placement.** Queries can attach to earlier-added queries rather than backbone edges. After pruning other queries the attachment edge is well-defined but reflects joint placement — document, don't "fix".
- **Backbone preservation [UNVERIFIED].** Wrapper test must assert output-tree-restricted-to-backbone has RF = 0 with the input tree.
- **Raw mode** (MinHash k=15) is unlikely to resolve 150 nt V4 fragments or SARS-CoV-2 genomes differing by a few SNPs. Primary DIPPER arm: `-i m` on the shared alignment; secondary: `-i r` as native end-to-end.
- **Maturity:** v0.1.x — pin exact version / image digest; expect CLI churn.

### 1.4 Transitive licensing

"Non-GPL" only matters if the lab **distributes** the tool (pipeline, QIIME 2 plugin, container, commercial product) or — for **AGPL** — **runs it behind a network service**. Running GPL tools internally on lab data imposes essentially no obligations.

| Candidate | Own license | Effective license of what executes |
|---|---|---|
| pplacer, SEPP, APPLES-2 | GPL-3 | GPL-3 (APPLES also calls FastTree, GPL-2+) |
| pplacer-SCAMPP | MIT | MIT wrapper + **GPL-3 engine** |
| EPA-ng-SCAMPP, BSCAMPP (default) | MIT | MIT wrapper + **AGPL-3 engine** — *stricter* than GPL for services |
| EPA-ng | AGPL-3 | AGPL-3 |
| DIPPER | MIT | MIT (Boost BSL, TBB Apache-2 [UNVERIFIED], CUDA runtime) |

Surrounding (non-placement) steps also carry licenses: RAxML-ng AGPL-3 [V], FastTree GPL-2+ [V], WITCH GPL-3 [V], HMMER BSD [V], MAFFT BSD [V]. **jplace → tree for UniFrac/Faith's PD** via `guppy tog` is GPL (pplacer); DIPPER emits a tree directly. For a fully non-copyleft *production path*, these steps matter as much as the placer.

---

## 2. Fairness / protocol decisions

| Decision | Recommendation | Challenge / why |
|---|---|---|
| **Query alignment** | **Two arms.** (A) *Placement-only*: every alignment-requiring tool and DIPPER `-i m` get **the same** query alignment, masked to backbone columns. (B) *End-to-end*: SEPP native, DIPPER `-i r`, and every other tool as `align_time + place_time` from A. | SEPP's runtime always includes alignment; comparing to pplacer-without-alignment flatters pplacer. The lab's real question ("replace SEPP") is end-to-end. |
| **Arm-A aligner** | RNAsim: **true alignment** (oracle). 16S: **WITCH**. SARS-CoV-2: reference-based to Wuhan-Hu-1 (`mafft --addfragments --keeplength`). | Single `hmmalign` profile is poor on diverse 16S backbones. WITCH is slow → time and report separately. Alternative: reuse SEPP's extended alignment (cheaper, but couples everything to SEPP). |
| **Backbone topology** | One topology per (dataset, size), shared by all tools. RNAsim: FastTree-2 on true alignment. 16S: **GG2 topology pruned** to subsample. SARS-CoV-2: MAT pruned. | RAxML-ng search impractical at 200k. |
| **Model parameters** | Topology fixed; parameters estimated **in each family's native format** on it: raxml-ng `--evaluate` GTR+G4 → EPA-ng family; FastTree GTR+CAT log → pplacer family; RAxML-8 `-f e` info → SEPP; FastTree `-nome -noml` → APPLES-2; DIPPER none. | "Same params for all" is impossible (incompatible formats, different branch-length semantics). Ref-prep is reported as one-time cost, excluded from placement time. |
| **Threads** | **16 physical cores/job** via SLURM `--cpus-per-task` + cgroup; every thread flag = 16. BSCAMPP `--threads 16 --cpus-per-job 2`. DIPPER: cgroup only. | Thread scaling out of scope. |
| **Query type** | **Full-length 16S** and **in-silico V4 fragments** cut from the same sequences; fragment length per Q4 (default **150 nt**, EMP Deblur). | Fragment length shifts accuracy more than tool choice → facet, don't pool. |
| **SEPP config** | Run SEPP **as the lab runs it in production** (Q2) *and* with defaults on shared backbones. | The bar is "as accurate as *our* SEPP". |
| **Placement reported** | Top-LWR placement per query. | Comparable across tools; DIPPER/APPLES only give one. |

---

## 3. Ground truth & accuracy metrics

### 3.1 Datasets and truth

| Dataset | Role | Truth | Caveat |
|---|---|---|---|
| **RNAsim** (1M simulated RNA; SCAMPP/BSCAMPP/APPLES benchmark; Dryad doi:10.5061/dryad.78nf7dq [link V, contents UNVERIFIED]) | **Primary accuracy** + scaling to 200k | True tree + true alignment | Simulated, not 16S. Reuse rather than simulate (KISS). |
| **GG2 full-length 16S backbone** (WoL2 genome-informed, uDance-built [V]) | Realistic 16S accuracy; full-length + V4 fragments | **Leave-one-out** vs. GG2 backbone position | Estimated tree → error floor. |
| **EMP / Qiita V4 ASVs** | **Runtime/memory at realistic scale** + **concordance with SEPP** | **None** (GG2 ASV positions were placed by DEPP [V]) | Report as "concordance", not accuracy. |
| **SARS-CoV-2** (UShER public MAT, MIT [V]) | Extreme high-similarity stress test | LOO on MAT | MAT *is* UShER placement → agreement-with-UShER. Polytomies/zero-length branches (§3.2). Different locus/length → separate panel (Q8). |

### 3.2 Metrics

- **Edge distance (primary):** edges between placed and true edge in the backbone. Collapse zero-length/unsupported branches first; placement on *any* edge of the true polytomy-equivalent set scores 0 (else SARS-CoV-2 is noise).
- **Delta error (secondary):** increase in FN bipartitions after insertion vs. reference restricted to backbone + query (Mirarab et al. 2012), for comparability with SCAMPP/BSCAMPP/APPLES papers.
- **Failures:** OOM / timeout / crash → **DNF**, plotted explicitly, never dropped.
- **"As accurate as SEPP" (Q9 default):** per bin with ≥100 queries, paired per-query difference (tool − SEPP) in edge distance; pass if upper 95% bootstrap CI of the mean ≤ 0.25 edges in every bin.

---

## 4. Defining "percent similarity"

**Per query:** `sim(q) = max_b identity(q, b)` over backbone leaves *b*, computed on the shared alignment over columns where *q* is non-gap (double gaps ignored). For fragments, identity over the V4 window. Dataset-level mean pairwise identity is reported only as a descriptor.

**Bins:** 75–80, 80–85, 85–90, 90–95, 95–97, 97–99, 99–99.9, ≥99.9. Below ~75% 16S alignment is unreliable → out of scope.

**Populating bins:**
1. **Natural sampling** — draw held-out queries, bin by `sim`. High bins fill naturally (SARS-CoV-2 fills ≥99%).
2. **Clade-pruning** for sparse low bins — for target [lo, hi), remove backbone leaves with identity ≥ hi to *q*, then **resubsample the backbone back to fixed size N** so bins share backbone size.

**Confounds to state in the report:**
- Pruning thins the local neighbourhood: low-similarity error measures "sparse neighbourhood", not just "divergent query".
- Low-identity queries over-represent certain lineages (archaea, CPR) → report per-bin lineage composition.
- **Similarity is confounded with backbone size** → always bin *within* a fixed backbone size, never pool across sizes.
- Fragments drift into higher bins than their parents → facet by query type.

---

## 5. Scale

| Axis | Values | Notes |
|---|---|---|
| Backbone leaves | **1k, 10k, 100k, 200k** | RNAsim covers all; GG2 up to its 16S count [UNVERIFIED]; SARS-CoV-2 1k/10k/100k |
| Queries (accuracy) | 1,000 per (dataset, size, query type), ≥100 per bin where possible | Under-filled bins flagged, not padded |
| Queries (scaling) | 1k, 10k, 100k on fixed 10k backbone | Matches ASV batch sizes (Q5) |
| Replicates | Accuracy: 3 query-set replicates. Timing: all runs timed; 3 repeats of one mid-size config for noise | Exclusive nodes, one node type |
| Limits | 256 GB / 48 h per job; exceed → DNF | |

**Expected walls** (all [UNVERIFIED] until T21 pilot): pplacer DNF at 100k–200k (SCAMPP's motivation); SEPP memory set by placement-subset size; EPA-ng control may DNF at 200k; DIPPER GPU bounded by GPU memory; APPLES-2, BSCAMPP, DIPPER should scale furthest.

---

## 6. Measurement

| Quantity | How | Gotcha |
|---|---|---|
| Wall-clock | Harness timestamps around tool call, cross-checked with `/usr/bin/time -v` | — |
| CPU time | `/usr/bin/time -v` user + sys (includes waited children) | — |
| **Peak memory** | **cgroup v2 `memory.peak`** (SLURM cgroup plugin or `systemd-run --scope`); `time -v` max RSS kept as secondary | `time -v` reports largest *single* process — under-reports forking tools. Snakemake `benchmark:` polls coarsely; not primary. |
| GPU memory | Poll `nvidia-smi --query-compute-apps=used_memory` every 0.5 s, record max | Approximate; record GPU model |
| Phase split | Separate rules: **ref-prep → align → place → convert → evaluate**, each measured | SEPP align+place inseparable → report total; SEPP rebuilds HMMs every run (real cost) |
| Provenance | CPU model, node, tool version, env lock hash, exact command | — |

Outputs — one tidy row per run (`results/runs.tsv`):
`dataset, backbone_n, query_type, query_n, bin, replicate, tool, arm, phase, wall_s, cpu_s, peak_mem_mb, gpu_peak_mb, status, threads, node, cpu_model, version`
and per-query accuracy (`results/accuracy.tsv`): `run_id, query, sim, edge_dist, delta_err`.

---

## 7. Stack (KISS)

| Concern | Choice | Why |
|---|---|---|
| Workflow | **Snakemake** (per-rule `conda:`, SLURM executor plugin) | Python-native; common in bio labs |
| Environments | **One pinned conda env YAML per tool** (bioconda/conda-forge, linux-64) + explicit lockfile; pip inside env for `apples`, `bscampp` | Conflicting deps isolated |
| Containers | Only **DIPPER GPU** (Apptainer `.sif` from Docker, pinned by digest) | HPC runs Apptainer |
| Analysis | Python 3.11, pandas, matplotlib, treeswift; hand-rolled jplace reader | No guppy/gappa in scoring |
| Tests | `pytest` for pure functions on hand-built trees; **smoke test**: Snakemake on committed tiny fixture (~100 leaves, <1 MB) running every tool (`dipper_cpu` in CI) on GitHub Actions `ubuntu-latest` (x86_64) | CI proves wiring; HPC proves scale |
| Lint | `ruff check`, `ruff format --check`, `snakemake --lint` | — |
| Compute | **HPC SLURM, x86_64**, exclusive nodes of one type; one GPU node type for DIPPER (also `dipper_cpu` on CPU nodes) | Q3 |
| Data | `data/`, `results/` gitignored; only fixture committed; scripted downloads with checksums | — |

---

## 8. Ticket breakdown & timeline

10 weeks, 5 milestones. See `TICKETS.md` for full ticket specs. Wrapper tickets T9–T15 parallelize once T8 is green; dataset tickets T16–T20 run in parallel with wrappers.

| Milestone | Weeks | Tickets |
|---|---|---|
| M0 Foundations | 1–2 | T1 skeleton/CI · T2 license & install audit · T3 fixture · T4 measurement harness · T5 accuracy metrics · T6 ref-prep · T7 query alignment · T8 tiny E2E smoke test |
| M1 Tool wrappers | 2–4 | T9 pplacer · T10 SEPP · T11 pplacer-SCAMPP · T12 EPA-ng-SCAMPP + EPA-ng · T13 BSCAMPP · T14 APPLES-2 · T15 DIPPER |
| M2 Datasets | 2–5 | T16 similarity/binning · T17 RNAsim · T18 GG2 16S · T19 EMP ASVs · T20 SARS-CoV-2 |
| M3 Runs | 5–8 | T21 HPC profile & pilot · T22 full benchmark |
| M4 Analysis | 8–10 | T23 plots · T24 report & recommendation |

Critical path: T1 → T2 → T6/T7 → T8 → T9–T15 → T21 → T22 → T23 → T24.

---

## 9. Open questions for the user

Each has a default used if unanswered.

1. **Why does GPL matter?** Internal use, redistribution in a pipeline/plugin, web service, or commercial? AGPL (EPA-ng → BSCAMPP, EPA-ng-SCAMPP) is *worse* than GPL for services. *Default:* internal use; license is a reported column, not a filter.
2. **Production SEPP setup?** q2-fragment-insertion? Which reference (GG 13_8, SILVA 128, GG2)? `-A`/`-P`? *Default:* q2-fragment-insertion defaults with its packaged reference, plus default SEPP on our backbones.
3. **Compute:** cluster, CPU model/cores, max memory/node, wall limit, GPU model/count? Is cgroup v2 `memory.peak` exposed? *Default:* SLURM x86, ≥16 cores, 256 GB, 48 h; DIPPER CPU-only if no NVIDIA GPU.
4. **Query length:** 150 nt V4, 250 nt, full-length, several? *Default:* 150 nt V4 + full-length.
5. **Typical unique ASVs per production placement run?** *Default:* 10k (scaling to 100k).
6. **16S backbone:** GG2, SILVA, GTDB? *Default:* GG2.
7. **Add plain EPA-ng control?** *Default:* yes.
8. **Keep SARS-CoV-2?** Different locus/length; UShER-circular truth. *Default:* yes, separate panel ≤100k, never pooled. Include UShER as reference point? *Default:* no.
9. **"As accurate as SEPP" threshold?** *Default:* per-bin paired non-inferiority margin 0.25 edges (§3.2).
10. **Deadline?** *Default:* 10 weeks from approval.
11. **Other candidates** (e.g. DEPP, which GG2 uses for V4 insertion [V])? *Default:* out of scope; revisit if nothing passes.
12. **Is GPU placement acceptable in production?** *Default:* rank on CPU numbers; GPU shown separately.
13. **Downstream needs inserted tree (UniFrac/Faith's PD)?** If yes, jplace→tree (guppy, GPL) is part of the replaced pipeline. *Default:* yes — report conversion cost and license.

---

## 10. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| pplacer/SEPP/EPA-ng hit memory/time walls at 100k–200k | Missing points | Expected; record DNF; T21 pilot calibrates |
| DIPPER: no jplace, joint placement, may alter backbone, v0.1.x churn, conda CPU-only | Bad/unfair scores | Adapter + RF check (T15); pin digest; report caveats |
| BSCAMPP uses bundled binaries / shared `~/.bscampp` | Wrong engine version; races | Per-job config + explicit binary paths (T2, T11–T13) |
| Circular truth (GG2 ASVs/DEPP; MAT/UShER) | Biased accuracy | Concordance/stress only; primary truth = RNAsim + GG2 LOO |
| Similarity confounded with density, clade, backbone size | Misread gradient | Bin within fixed size; resubsample after pruning; per-bin lineage composition |
| Memory under-reported for forking tools | Wrong conclusions | cgroup `memory.peak` primary |
| Timing noise / node heterogeneity | Spurious differences | Exclusive single-type nodes; repeats; record CPU model |
| Alignment quality dominates fragment accuracy | Tool effects masked | Oracle-alignment arm on RNAsim; shared alignment in arm A |
| Expensive ref-prep (RAxML-8 info at 200k for SEPP) | Schedule slip | Start T6 early on largest backbones; one-time cost |
| Fresh upstream releases (pplacer alpha23, epa-ng 0.3.9-pre) | Irreproducibility | Pin + lock every env |
| Storage (RNAsim 1M, MAT, EMP ASVs) | Disk pressure | Subsample at ingest; gitignored `data/`; checksummed downloads |
| Shared aligner (WITCH) is GPL and slow | Inflated end-to-end cost | Report alignment time separately |

---

## Handoff notes for SWE adversarial review

Attack in particular: (a) are arms A+B the minimum design that answers "replace SEPP"? (b) is cgroup `memory.peak` available on the target cluster? (c) WITCH vs. SEPP's extended alignment as the shared alignment; (d) is 24 tickets too granular (merge T3/T5 or T16/T17)? (e) is DIPPER's prune-others scoring fair to joint placement?
