# Simulated 10k true-tree panel (T16, wired in T22)

A uDance-free accuracy cross-check: GG2 truth was seeded by uDance (APPLES-2), so apples2 results on GG2 carry a caveat (`design.md` section 4). Here the truth is the simulated tree itself.

## Build

`sim10k_panel` (`workflow/rules/sim10k.smk`) runs `workflow/scripts/make_sim_panel.py` (T3 simulator, fixed seed) into `data/sim10k/B10k/` and moves the query files up one level, where the shared rules expect `data/{dataset}/queries*.fasta`:

| File | Content |
|---|---|
| `data/sim10k/B10k/{tree.nwk, aligned-dna-sequences.fasta, raxml-info.txt}` | backbone |
| `data/sim10k/B10k/reference_tree.nwk` | truth: the full simulated tree |
| `data/sim10k/B10k/truth.tsv` | held-out clades and singletons |
| `data/sim10k/queries.fasta`, `queries.aligned.fasta` | held-out queries (the aligned file is the oracle alignment, unused here) |

The panel is 10,000 simulated leaves minus 1,025 held-out queries (7 x (5+10+15+20+25) clade members and 500 singletons), so the backbone has **8,975 leaves**. It is named `B10k` so the shared rules and the per-backbone binning treat it like the GG2 B10k\*, but it is not the same tree: never pool it with `gg2/B10k`.

## Cells

Backbone `B10k` only (no Bfull\*), regions `full-length` and `V4-150`, `hmmalign` only (no oracle), all 11 arms: 22 cells (`exclude` in `config/full.yaml`). The simulator embeds primer sites, so every query survives V4-150 extraction (`tests/test_sim_panel.py`).
