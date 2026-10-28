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

`V4-150c`, `V3-V4c`, `V1-V2c`, `V1-V3c`, `V4-V5c`, `V5-V6c`, `V6-V8c`, `V7-V9c` cut a query by alignment columns instead of primer sites: the **column-window proxy** (real amplicons start and stop at primers, so these are not amplicons). The primer-matched regions are unchanged and never pooled with these.

- Window: over the shipped aligned query rows that have both primer hits, the median start and median end column of the primer-delimited amplicon (`regions.column_window`; V4-150c is the first 150 nt after 515F expressed in columns). Windows come from the queries, not a backbone, because region query files are dataset-level.
- Cut: the query's own shipped aligned row, columns `[start, end)`, gaps removed (`regions.cut_column_window`); `results/gg2/regions/<region>/queries.fasta` then goes through the same hmmalign, place and evaluate rules as any region.
- Kept when its non-gap residues in the window are at least 0.9 x the window length (median residue count over the primer-hit rows); the fraction is `COVERAGE_TOLERANCE`, fixed before the run and the same for every region. The rest are listed in `exclusions.tsv` and counted in `exclusion_counts.tsv` as `window coverage below tolerance`.
- A query without primers that covers the window is kept. That is the point: 711 of the 1,000 sources are MJ0* UMI amplicons, and V1-V2, V1-V3 and V7-V9 keep only 153, 144 and 46 queries by primers. Why those reads lack the primer sites (the UMI amplicons being primer-trimmed) is UNVERIFIED pending the `exclusion_counts.tsv` check in `docs/runbook-t28.md` stage S0.
- A window is a position-based proxy: a query that starts or ends inside the window by a different amount than a real amplicon would is still cut at the same columns, so accuracy is comparable across arms, not to amplicon studies.
- Not planned: the oracle alignment (it is defined from primer spans), and any EMP cell.

## Licenses

GG2 data license: BSD-3 [UNVERIFIED], flagged only (`LICENSES.md` is unchanged).

## Running on the HPC

```sh
srun -c 16 --mem=32G snakemake --use-conda -c 16 data/gg2/MANIFEST.tsv
```

Not run on macOS (about 8.5 GB of similarity memory and `hmmalign` are not run locally).
