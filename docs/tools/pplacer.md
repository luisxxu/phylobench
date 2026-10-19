# pplacer (`place`, arm `pplacer`)

Baseline placement arm. Rule: `workflow/rules/pplacer.smk` (`pplacer_place`); arm file `config/arms/pplacer.yaml`. This is the template for new wrappers (`docs/harness.md` section 7).

## Command line

```sh
phylobench-measure --out {cell}/place.run.json --threads {threads} --timeout 86400 \
    --tool-version 1.1.alpha23 --env-lock envs/pplacer.linux-64.pin.txt -- \
    pplacer -j {threads} -t data/{dataset}/{backbone}/tree.nwk \
    -s data/{dataset}/{backbone}/raxml-info.txt \
    -r data/{dataset}/{backbone}/aligned-dna-sequences.fasta \
    -o {cell}/placements.jplace results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta
```

stdout/stderr go to `{cell}/place.log`. On `status != ok` the rule touches `placements.jplace` so `run.json` survives (`docs/harness.md` section 2). The `convert` phase (`guppy tog`) follows in the same cell (`docs/tools/guppy.md`).

`--timeout 86400` (24 h) applies to every pplacer cell, a deviation from `design.md:124` (which caps only the Bfull\* attempt); a timeout is a `status=timeout` DNF row. A timed-out cell also gets a `convert` DNF error row, because `guppy tog` runs on the empty jplace.

## Version

Pinned `pplacer=1.1.alpha23`, build `h01ede22_0` (`envs/pplacer.yaml`, lock `envs/pplacer.linux-64.pin.txt`). `pplacer --version` prints only `dev`; the proof is the package metadata (`docs/audit.md` section 2). `--tool-version 1.1.alpha23` is hard-coded in the rule; update it with the pin.

## Threads

`-j {threads}` (config `threads`, 16 in production) plus the SLURM reservation and the harness `threads` field.

## License

GPL-3, and the engine executed is pplacer itself (see `LICENSES.md`).

## Capped Bfull\* attempt

Design section 5 plans one capped pplacer attempt at Bfull\* (1k V4 queries, 24 h cap, no retry). Run it with the overlay:

```sh
snakemake --use-conda --profile profiles/slurm --configfile config/full.yaml config/pplacer-bfull-cap.yaml
```

`config/pplacer-bfull-cap.yaml` names a placeholder dataset, `gg2-v4-1k`, that T16/T21 data prep must produce. There is no retry: the harness exits 0 on a DNF, and T21's profile must keep retries at 0. `runs.tsv` is regenerated on every invocation, so the final full run regenerates it.

## Provenance

The pplacer `url#md5` line in `envs/pplacer.linux-64.pin.txt` is identical across the pplacer, sepp and bscampp locks (the lock-guard in `.github/workflows/envs.yml` enforces it), so all pplacer-based arms share one binary. Every `run.json` records `env_lock_sha256`.

## Quirks

- alpha20 deadlocks with `-j`; alpha23 is the pinned fix.
- `-s` reads the RAxML-8 info file shipped with the reference (model numbers are the shipped ones, never re-estimated).
- Aligned arm: it places the `hmmalign` or `oracle` aligned queries, which have the backbone's column count (checked by `check_smoke.py`).
