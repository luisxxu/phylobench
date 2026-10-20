# EMP V4 ASV query sets (T19)

Real Earth Microbiome Project Deblur V4 ASVs, used for runtime/memory scaling and for concordance with `sepp-prod`.

**EMP positions are not ground truth (DEPP circularity: GG2 placed these ASVs with DEPP). Outputs are labelled `concordance`, never `accuracy`.**

## Source

| | |
|---|---|
| URL | `https://ftp.microbio.me/emp/release1/otu_tables/deblur/emp_deblur_150bp.release1.biom` |
| Size | 167,259,648 bytes |
| Last-Modified | 2017-09-12 |
| sha256 | `b324e5098666839ab9747bd7d78e5c7a2f81cc3f249eaeea5ccd755b6c18196a` |
| Format | HDF5, BIOM format-version `[2, 1]` |
| Shape | 197,711 observations x 17,483 samples |
| Observation IDs | 197,711, all unique, all match `^[ACGT]{150}$` |

No checksum is published for this file. The sha256 above was pinned on first download and is enforced by `workflow/scripts/fetch.sh` (shared with the GG2 download, `docs/gg2.md`) on every later download.

## Query sets

- `emp-10k` and `emp-100k`: uniform samples of the unique ASVs, seed 0 (`phylobench.emp.query_set`: sorted IDs, `numpy.random.default_rng(0).permutation`, first n).
- The sets are nested: `emp-10k` is the first 10,000 ASVs of `emp-100k`.
- ASVs are renamed `asv000001`, `asv000002`, ... because 150-character names are fragile across tools. The name-to-ASV mapping is the FASTA order of `results/{dataset}/regions/V4-150/queries.fasta`. The script fails if a name collides with a backbone leaf.
- ASVs are already V4 amplicons without primers, so the FASTA is written directly as region `V4-150` (no primer matching). Similarity is computed over each query's own non-gap columns, so no window is needed.
- `alns: [hmmalign]` only: EMP has no oracle alignment.

## Concordance with sepp-prod

`concordance` rule: for each EMP cell, `edge_dist` between the arm's placement and `sepp-prod`'s placement of the same ASV on the same backbone (T5 edge distance with the `sepp-prod` node as the reference position). Output `results/concordance.tsv`, `label=concordance`. Neither placement is truth. Bins are never pooled across backbones.

`sepp-prod` is also an arm of the EMP config: its own rows score against themselves (distance 0) and must be filtered out by consumers (`run_id == reference_run_id`).

## Pipeline

- Run the gg2 stage first. EMP reuses GG2's backbones (`Bfull`, `B10k`): `data/emp-*/{backbone}/{tree.nwk,aligned-dna-sequences.fasta,raxml-info.txt}`, `results/ref-prep/emp-*/{backbone}/{me-tree.nwk,epa-ng.model}` and `results/emp-*/{backbone}/align/backbone.hmm` are symlinks to the gg2 files. There is one ref-prep and one HMM per backbone.
- `data/emp/MANIFEST.tsv` (`dataset kind backbone region bin path n sha256 flagged`): the biom, each queries FASTA, each linked backbone file (sha256 verified against `data/gg2/MANIFEST.tsv`), and one `kind=bin` row per dataset x backbone x bin, zero-count bins included, `flagged=1` when n < 100.
- Config: `config/emp.yaml`. `data/emp/MANIFEST.tsv` is only requested when a configured dataset is an EMP dataset, so the smoke test never downloads.

## GG2 interface used (confirmed against T18)

- `workflow/scripts/fetch.sh <url> <sha256> <out>`, called directly through `params.fetch` as in `gg2_download`.
- Backbone names `Bfull`, `B10k`; `data/gg2/{Bfull,B10k}/{tree.nwk,aligned-dna-sequences.fasta,raxml-info.txt}` come from `gg2_build`, so requesting them plans `gg2_unzip` and `gg2_build` (never `gg2_download` once the `.qza` files exist).
- `data/gg2/MANIFEST.tsv` columns `kind backbone region bin path n sha256 flagged`; `kind=file` rows are keyed by repo-relative `path`, and `raxml-info.txt` has an empty `n`. `emp_manifest.py` reuses `phylobench.gg2.file_row` and `bin_rows` and adds a `dataset` column.
- gg2 rules use literal `data/gg2/...` outputs, so they never match the `emp-10k|emp-100k` wildcard rules.
- Standard ref-prep and HMM paths: `results/ref-prep/gg2/{backbone}/{me-tree.nwk,epa-ng.model}`, `results/gg2/{backbone}/align/backbone.hmm`.
