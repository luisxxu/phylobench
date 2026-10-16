# pplacer (`place`, arm `pplacer`)

Baseline placement arm. Rule: `workflow/rules/pplacer.smk` (`pplacer_place`); arm file `config/arms/pplacer.yaml`. This is the template for new wrappers (`docs/harness.md` section 7).

## Command line

```sh
phylobench-measure --out {cell}/place.run.json --threads {threads} \
    --tool-version 1.1.alpha23 --env-lock envs/pplacer.linux-64.pin.txt -- \
    pplacer -j {threads} -t data/{dataset}/{backbone}/tree.nwk \
    -s data/{dataset}/{backbone}/raxml-info.txt \
    -r data/{dataset}/{backbone}/aligned-dna-sequences.fasta \
    -o {cell}/placements.jplace results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta
```

stdout/stderr go to `{cell}/place.log`. On `status != ok` the rule touches `placements.jplace` so `run.json` survives (`docs/harness.md` section 2). The `convert` phase (`guppy tog`) follows in the same cell (`docs/tools/guppy.md`).

## Version

Pinned `pplacer=1.1.alpha23`, build `h01ede22_0` (`envs/pplacer.yaml`, lock `envs/pplacer.linux-64.pin.txt`). `pplacer --version` prints only `dev`; the proof is the package metadata (`docs/audit.md` section 2). `--tool-version 1.1.alpha23` is hard-coded in the rule; update it with the pin.

## Threads

`-j {threads}` (config `threads`, 16 in production) plus the SLURM reservation and the harness `threads` field.

## License

GPL-3, and the engine executed is pplacer itself (see `LICENSES.md`).

## Quirks

- alpha20 deadlocks with `-j`; alpha23 is the pinned fix.
- `-s` reads the RAxML-8 info file shipped with the reference (model numbers are the shipped ones, never re-estimated).
- Aligned arm: it places the `hmmalign` or `oracle` aligned queries, which have the backbone's column count (checked by `check_smoke.py`).
