# SEPP production (`align+place+convert`, arm `sepp-prod`)

Yardstick arm: production SEPP settings and versions, run standalone outside QIIME 2. Rule: `workflow/rules/sepp-prod.smk` (`sepp_prod`); arm file `config/arms/sepp-prod.yaml`; converter `workflow/scripts/sepp_jplace.py` (`phylobench.jplace.from_sepp`).

## Command line

Run from a per-job work dir `{cell}/work` with `TMPDIR={cell}/work/tmp`:

```sh
phylobench-measure --out {cell}/sepp.run.json --threads {threads} \
    --tool-version 'sepp 4.5.6, pplacer 1.1.alpha19' --env-lock envs/sepp-prod.linux-64.pin.txt -- \
    run-sepp.sh results/{dataset}/regions/{region}/queries.fasta sepp -x {threads} -A 1000 -P 5000 \
    -a data/{dataset}/{backbone}/aligned-dna-sequences.fasta \
    -t data/{dataset}/{backbone}/tree.nwk -r data/{dataset}/{backbone}/raxml-info.txt
```

The cell is `results/{dataset}/{backbone}/{region}/raw/sepp-prod/`. stdout/stderr and SEPP's own `sepp-sepp-*.log` go to `sepp.log`. On success `work/sepp_placement.json` is converted to `placements.jplace` (unmeasured, `from_sepp`); on `status != ok` the rule writes an empty `placements.jplace` so `run.json` survives.

## Version

Own env `envs/sepp-prod.yaml` and lock `envs/sepp-prod.linux-64.pin.txt`: sepp 4.5.6, pplacer 1.1.alpha19, hmmer 3.4, openjdk 25.0.2, python 3.12.13: the versions QIIME 2 2026.7 ships (`envs/qiime2.linux-64.pin.txt`). Not QIIME 2 itself, and not in the `envs.yml` install matrix. Builds differ from the QIIME 2 lock where the yaml cannot pin builds (the lock-guard regex forbids it):

| Package | `sepp-prod` lock | QIIME 2 2026.7 lock |
|---|---|---|
| sepp | `4.5.6-py312h87e0c26_4` | same |
| pplacer | `1.1.alpha19-h9ee0642_2` | same |
| hmmer | `3.4-h7d74f8d_5` | same |
| openjdk | `25.0.2-h1602c4f_26` | `25.0.2-ha668962_0` |
| python | `3.12.13-h8ab3286_1_cpython` | `3.12.13-hd63d673_0_cpython` |

## Threads

`-x {threads}` (config `threads`, 16 in production) plus the SLURM reservation and the harness `threads` field.

## Work and temp dirs

One work dir per job (`{cell}/work`, wiped at start, `.tog` files kept) and `TMPDIR` set to `{cell}/work/tmp`. `run-sepp.sh` overrides `-d` and `-p` itself (a `mktemp` dir), so `TMPDIR` is what isolates temps between jobs.

## License

SEPP GPL-3, pplacer GPL-3 (including `guppy`), HMMER BSD-3, OpenJDK (GPLv2 with Classpath exception). The engines executed are SEPP, pplacer and HMMER (see `LICENSES.md`).

## Quirks

- SEPP writes jplace version 1 with `:len[N]` edge numbers; `from_sepp` rewrites them to version 3 `{N}` and fails if there is no number, a `{` already present, or a missing required field. Internal nodes are renamed by SEPP; scoring maps edges by splits.
- SEPP resolves polytomies in the backbone, so the GG2 backbone must be binary (T18).
- Fragments SEPP cannot score are dropped with only a log line; they surface as unplaced queries in `evaluate`.
- SEPP's default seed 297834 (`sepp/config.py`); `run-sepp.sh` doesn't pass `-seed`.
- `guppy tog` runs inside the phase, as in production; the single phase covers alignment, placement and conversion.
- The benchmarked cell uses the pruned backbone (Bfull\*/B10k\*), not the as-shipped reference.
