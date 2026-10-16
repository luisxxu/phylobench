# HMMER (hmmbuild + hmmalign) — shared query alignment

Pinned `hmmer=3.4` (`envs/hmmer.yaml`, lock `envs/hmmer.linux-64.pin.txt`). Rules: `workflow/rules/align.smk`. Both calls run through `phylobench-measure` in phase `align` (`results/{dataset}/{backbone}/align/hmmbuild.run.json`, `.../{region}/align/hmmalign.run.json`), with the DNF contract from `docs/harness.md` §2.

## Pipeline per backbone and region

1. `backbone_stockholm`: backbone alignment FASTA to one-block Stockholm with an all-`x` `#=GC RF` line.
2. `hmmbuild --dna --hand --cpu {threads} backbone.hmm backbone.sto`. `--hand` takes the match columns from the RF line, so every backbone column becomes a match state and the profile has exactly the backbone's column count. Without `--hand` hmmbuild drops gappy columns and the masked width would not match.
3. `extract_regions`: cut the ungapped full-length queries by primer matching (`src/phylobench/regions.py`); also writes `exclusions.tsv` and `exclusion_counts.tsv` per region.
4. `hmmalign --dna --outformat afa -o hmmalign.afa backbone.hmm queries.fasta` (queries only; no `--mapali`).
5. `mask_hmmalign`: delete insert columns (lowercase residues and `.` gaps) with `str.translate`, giving `hmmalign.aligned.fasta`. The script fails loudly if any masked row is not the backbone's width. T8 asserts the same end to end.

`oracle_region` (B10k\* diagnostic) writes `oracle.aligned.fasta`: the query's own GG2 alignment columns for the region, gaps elsewhere, same width as the backbone. It first checks that the aligned rows ungap to the full-length FASTA.

## Inputs (not committed)

`data/{dataset}/queries.fasta` (unaligned full-length), `data/{dataset}/queries.aligned.fasta` (GG2 rows, oracle only), `data/{dataset}/{backbone}/aligned-dna-sequences.fasta` (pruned backbone alignment).

## Threads

`hmmbuild` gets `--cpu 16` (config `threads`). `hmmalign` 3.4 has no `--cpu` option (checked with `hmmalign -h`), so it is single-threaded; the rule still requests 16 threads from SLURM and records `threads` in `run.json`, so the scheduler allocation matches every other arm.

## Regions

Primers are matched exactly (IUPAC-aware, 0 mismatches) on the ungapped sequence; the amplicon excludes both primers. Sequences with no match are excluded and counted per reason (`no forward primer`, `no reverse primer`, `too short`).

| Region | Forward | Reverse |
|---|---|---|
| V4-150 | 515F (Parada) `GTGYCAGCMGCCGCGGTAA` | 806R (Apprill) `GGACTACNVGGGTWTCTAAT` |
| V3-V4 | 341F `CCTACGGGNGGCWGCAG` | 805R `GACTACHVGGGTATCTAATCC` |
| V1-V2 | 27F `AGAGTTTGATCMTGGCTCAG` | 338R `GCTGCCTCCCGTAGGAGT` |

V4-150 is the first 150 nt after 515F, and 806R must still match downstream. A DNF'd `hmmbuild` leaves an empty `.hmm`, so the dependent `hmmalign` fails loudly rather than producing a row; a DNF'd `hmmalign` yields an empty masked output and its `run.json` carries the status.
