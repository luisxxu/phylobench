# GG2 2022.10 dataset (T18)

Builds the shared inputs of every GG2 run: Bfull\*, B10k\*, 1,000 held-out queries, per-region similarity bins. Rules: `workflow/rules/gg2.smk` (`gg2_download`, `gg2_unzip`, `gg2_build`); pure logic: `src/phylobench/gg2.py`, `src/phylobench/heldout.py`. All outputs are literal `data/gg2/...` paths (not a `{dataset}` wildcard) and are gitignored.

## Provenance

| File | URL | sha256 (pinned in `gg2.smk`) |
|---|---|---|
| `2022.10.backbone.sepp-reference.qza` (77,441,606 B) | `https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.sepp-reference.qza` | `fa5d32f936a39735a4f99552044ed4fcbae0d91fb48a679f9d329e2a708eb973` |
| `2022.10.backbone.tax.qza` (4,721,150 B) | `https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.tax.qza` | `6495d6da4097ff51391ae7ce8df09c0ee00807c8bb64e4ce8cae1c9d5e07f5f0` |

Both hashes were measured on the first download (2026-10-20) and are pinned; `workflow/scripts/fetch.sh <url> <sha256> <out>` removes the download and fails on a mismatch. The `raxml-info.txt` inside the reference archive has sha256 `1ad6d7b1b2c696a93ff5dcf2be1e677817bb50063474e41c0c0d3dfc22449f76`.

Members read (`<uuid>/data/`): reference `tree.nwk`, `aligned-dna-sequences.fasta`, `raxml-info.txt`; taxonomy `taxonomy.tsv`. `workflow/scripts/unzip_qza.py` extracts only these (stdlib `zipfile`, no QIIME 2).

## Measured properties of the shipped reference

- R = 331,269 leaves and alignment rows; alignment width L = 1,598 columns.
- Tree: strictly bifurcating (2-child root, every other internal node has 2 children); leaf labels equal the alignment names.
- Exact duplicate aligned rows (upper-case, `.` as `-`): 0 in the shipped alignment. `gg2.check_backbone` still enforces none in either backbone (FastTree needs it, `docs/tools/fasttree.md`).
- `taxonomy.tsv`: tab-separated, header `Feature ID<TAB>Taxon`, ranks `d__...; p__...; ...; s__...` joined by `; `, empty ranks keep their prefix (`g__`). Its IDs spell the tree's underscores as dashes (`RS-GCF-...` vs `RS_GCF_...`); the mapping is one to one, and `gg2.parse_lineages` maps them back.
- Similarity memory is about 16 * R * L bytes = 8.47 GB (8,078 MiB) for the full backbone; `gg2_build` requests 32,000 MB and `build_gg2.py` fails fast if 16 * R * L exceeds `--mem-mb`.

## Procedure

1. Held-out groups from the shipped tree and `taxonomy.tsv` (`heldout.choose_groups`, seed 2022): disjoint monophyletic clades with 25 to 2,000 leaves, taxa per rank (phylum 2, class 3, order 4, family 6, genus 6), then singletons. Clade queries are up to 25 random members per clade; the singleton count is 1,000 minus the clade queries (525 clade queries and 475 singletons with these settings).
2. Bfull\* = shipped minus every held-out leaf (`heldout.prune_tree`, `prune_alignment`), then exact duplicate rows removed (D1).
3. B10k\* = seeded uniform subsample of Bfull\* of 10,000 leaves, all alignment columns kept.
4. `raxml-info.txt` is a byte copy of the shipped file in each backbone dir (never re-estimated).
5. `reference_tree.nwk` per backbone = the shipped tree restricted to that backbone plus the 1,000 queries; `accuracy.true_positions` reads queries' true edges from it.
6. Queries: `queries.fasta` (ungapped) and `queries.aligned.fasta` (the shipped rows, full width). Region windows come from `regions.amplicon_span` on the ungapped query projected onto its aligned row; a query lacking a region's primers is not binned for that region.

## >=99.9 seeding (D2)

Random singletons in a 331k backbone rarely have a >=99.9 neighbour. Seeding keeps the clades and the random singleton stream identical and replaces only the trailing singletons:

1. Draw the groups; count singletons whose full-length similarity to the remaining backbone is >=99.9 (`hits`).
2. If the count is at least 100, stop.
3. Otherwise draw 2,000 candidate leaves of the remaining backbone (own RNG, `default_rng([2022, 1])`); candidates whose similarity to the backbone minus all candidates is >=99.9 are the seeds, in draw order.
4. `k = seeds_needed(hits, ...)`: the smallest k with `sum(hits[:n-k]) + k >= 100`; fail if k exceeds the seeds found (raise `N_CANDIDATES`).
5. Redraw with the first k eligible seeds taken first (no RNG), then the remaining random singletons.

A seed's neighbour can itself be held out by the redraw, so the full-length >=99.9 bin count of each built backbone is logged; fewer than 100 is flagged in `MANIFEST.tsv` (and logged), not an error. `heldout.tsv` marks seeded leaves `seeded=1`.

## Outputs

`data/gg2/{Bfull,B10k}/{tree.nwk,aligned-dna-sequences.fasta,raxml-info.txt,reference_tree.nwk}`, `queries.fasta`, `queries.aligned.fasta`, `heldout.tsv`, `MANIFEST.tsv`. Trees are written with `simulate.format_newick` (leaf labels only, 6-decimal lengths, no root length).

`heldout.tsv` columns: `leaf group group_type rank taxon query seeded`; one row for every held-out leaf, including clade members that are not queries (`query` 0 or 1).

`MANIFEST.tsv` columns, in order: `kind backbone region bin path n sha256 flagged`.

- `kind=file`: `path` (repo-relative), `n` (leaves or records where meaningful), `sha256`; `backbone` set for per-backbone files; `region`, `bin`, `flagged` empty. One row per output plus both `.qza` and `heldout.tsv`.
- `kind=bin`: `backbone` in {Bfull, B10k}, `region` one of the nine regions, `bin` a `similarity_bin` label (empty bins included, `n=0`), `n`, `flagged=1` when `n < 100`; `path` and `sha256` empty.

Bins are never pooled across backbones. These counts are advisory (window similarity against the built backbone); `results/similarity.tsv` from the workflow is authoritative. At B10k\* a small >=99.9 bin is expected (sparser neighbourhoods) and is flagged, not padded.

## Deviations from the first plan

- D1: exact duplicate aligned rows are removed from Bfull\* after hold-out (FastTree needs none); recorded in the log.
- D2: >=99.9 seeding replaces trailing singletons with candidates instead of adding extra ones, so the query count stays 1,000.
- D6 (polytomies in the shipped tree) did not apply: the shipped tree is strictly bifurcating. D3: manifest bin counts come from the oracle (shipped) alignment and are advisory; `results/gg2/{backbone}/{region}/similarity.tsv` (hmmalign) is authoritative. D4: both qza sha256 values are pinned from a one-time download (no upstream checksums are published). D5: `accuracy.true_positions` accepts extra non-backbone, non-query leaves in the reference tree (needed when primer exclusions drop queries from a region).

## Column-window regions (T28)

`V4-150c`, `V3-V4c`, `V1-V2c`, `V1-V3c`, `V4-V5c`, `V5-V6c`, `V6-V8c`, `V7-V9c` cut a query by alignment columns instead of primer sites: the **column-window proxy** (real amplicons start and stop at primers, so these are not amplicons). The primer-matched regions are unchanged and never pooled with these. They apply to `gg2` and `gg2-wol` (the rule reads `data/<dataset>/queries.aligned.fasta`); T28 uses only `V7-V9c` on `gg2-wol`, where the primer-matched V7-V9 keeps about 139 of 3,000 queries.

- Window: over the shipped aligned query rows that have both primer hits, the median start and median end column of the primer-delimited amplicon (`regions.column_window`; V4-150c is the first 150 nt after 515F expressed in columns). Windows come from the queries, not a backbone, because region query files are dataset-level.
- Cut: the query's own shipped aligned row, columns `[start, end)`, gaps removed (`regions.cut_column_window`); `results/gg2/regions/<region>/queries.fasta` then goes through the same hmmalign, place and evaluate rules as any region.
- Kept when its non-gap residues in the window are at least 0.9 x the window length (median residue count over the primer-hit rows); the fraction is `COVERAGE_TOLERANCE`, fixed before the run and the same for every region. The rest are listed in `exclusions.tsv` and counted in `exclusion_counts.tsv` as `window coverage below tolerance`.
- A query without primers that covers the window is kept. That is the point: 711 of the 1,000 sources are MJ0* UMI amplicons, and V1-V2, V1-V3 and V7-V9 keep only 153, 144 and 46 queries by primers. Why those reads lack the primer sites (the UMI amplicons being primer-trimmed) is UNVERIFIED pending the `exclusion_counts.tsv` check in `docs/runbook-t28.md` stage S0.
- A window is a position-based proxy: a query that starts or ends inside the window by a different amount than a real amplicon would is still cut at the same columns, so accuracy is comparable across arms, not to amplicon studies.
- Not planned: the oracle alignment (it is defined from primer spans), and any EMP cell.

## WoL genome held-out set (T31)

`gg2-wol` is a second held-out set on the same shipped reference, for the T28 to T30 follow-up: queries come only from leaves labelled `^G\d{9}$` (WoL genome leaves, `gg2.WOL_LEAF`). Rules: `gg2_wol_build` in `workflow/rules/gg2.smk` (literal `data/gg2-wol/...` paths, own `MANIFEST.tsv`, reuses the `gg2_unzip` outputs); `build_gg2.py --wol`; pure logic in `gg2.py` (`wol_options`, `wol_n_queries`) and `heldout.py` (`eligible`, `min_eligible`). The dataset key is `gg2-wol`; backbones keep the names `Bfull` and `B10k`, and it is never pooled with `gg2` (T22 stays as published). The T23/T24 figures and report filter `dataset == "gg2"`, so `gg2-wol` is reported by the T28 and T30 close-outs.

- Truth wording: **G-leaf genome position in GG2** (the leaf's place in the shipped tree). Whether uDance kept the WoL2 topology fixed when it built GG2 is UNVERIFIED; the apples2 caveat of design.md section 3 stays.
- Draw (`choose_seeded_groups(..., eligible=G leaves, min_eligible=5)`, seed 2026): disjoint monophyletic clades of 25 to 2,000 leaves with at least 5 G members, held out whole (G and non-G members both removed), at most 25 G members per clade as queries; clade counts per rank are `WOL_RANKS` (phylum 6, class 9, order 12, family 18, genus 18), defaults that the S0 probe (P9) finalizes. Singletons and the >=99.9 seeds are G leaves only; the singleton count is N minus the clade queries, counted as `min(G members, 25)` per clade. Filtering only: no extra RNG draws, so `eligible=None` reproduces the T22 draw exactly (regression test in `tests/test_heldout.py`).
- N = 3,000 queries; 2,000 if there are 6,000 to 9,999 G leaves; below 6,000 `wol_options` raises (G1 in `docs/runbook-t31.md`).
- Bfull\* = shipped minus every held-out leaf, duplicate rows removed; B10k\* = seeded uniform 10,000-leaf subsample of that Bfull\* (seed 2026). Per-backbone ref-prep (FastTree me-tree, epa-ng model, hmmbuild) runs fresh for `gg2-wol`.
- Outputs: `data/gg2-wol/{Bfull,B10k}/{tree.nwk,aligned-dna-sequences.fasta,raxml-info.txt,reference_tree.nwk}`, `queries.fasta`, `queries.aligned.fasta`, `heldout.tsv`, `MANIFEST.tsv`; same schemas as `data/gg2`. Not in `config/smoke.yaml`.
- S0 probe: `workflow/scripts/probe_wol.py` (helpers in `src/phylobench/wol.py`) is read-only and prints P1 to P10 and the G1 to G5 verdicts; the runbook has the pre-registered stop/go and every cluster command.

### Build (2026-11-09, `db699d2`)

- 3,000 queries, all `^G\d{9}$`: 2,195 singletons and 805 from 63 clades (phylum 6, class 9, order 12, family 18, genus 18); 0 seeded (the Bfull\* full-length >=99.9 bin already held 234 queries).
- Primer-matched queries per region (both backbones): full-length 3,000, V4-150 2,703, V3-V4 2,542, V4-V5 2,622, V6-V8 2,354, V1-V2 2,119, V1-V3 2,058, V5-V6 1,913, V7-V9 139.
- Full-length bins, Bfull\*: <75 20, 75-80 134, 80-85 169, 85-90 283, 90-95 482, 95-97 407, 97-99 631, 99-99.9 640, >=99.9 234. B10k\*: 67, 267, 322, 489, 864, 415, 382, 165, 29. V4-150 Bfull\*: 6, 54, 158, 158, 322, 160, 249, 184, 1,412. B10k\*: 41, 173, 313, 328, 596, 312, 406, 166, 368. Bins under 100 are flagged, not padded.
- Ref-prep and hmmbuild for both backbones: 8 jobs, all `ok`, 17.7 min wall, 0.1 node-hours by `sacct`.
- sha256: `MANIFEST.tsv` `cfd769508e8699376391fac4cfd2ac1621286c2dd32d0c964fe25107758f0dbf`, `heldout.tsv` `1709493399fd86eee875dc18ec4d0f76f8118465e038db77600d89a26bac7543`, `queries.fasta` `da8eb9b2ab2c9117e687454eb566478cdd378621acdfa9eb64fefa84eca512e3`, `Bfull/tree.nwk` `f04080d060b68e0369557f0dfcd1bb867bb4ff71ff5296cb64a5de0d47e03ce8`, `B10k/tree.nwk` `31f3c13af36139ec56d0ed0551423b7629eb3f134e0f782f416ab87b65ae9160`.
- Archive: `/ddn_scratch/lxxu/phylobench-archive/t31-gg2-wol.tgz` (77 MB; `data/gg2-wol`, T31 logs), sha256 `c8251fe3e68c41a685c93ca60d289b0ee80184307cb2a7e40e5b26a0a477dc91`.

## Licenses

GG2 data license: BSD-3 [UNVERIFIED], flagged only (`LICENSES.md` is unchanged).

## Running on the HPC

```sh
srun -c 16 --mem=32G snakemake --use-conda -c 16 data/gg2/MANIFEST.tsv
```

Not run on macOS (about 8.5 GB of similarity memory and `hmmalign` are not run locally).
