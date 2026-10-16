# guppy (`guppy tog`, phase `convert`)

Converts a jplace to an inserted tree for every jplace arm (`sepp-prod` runs it internally, as production does). Rule: `workflow/rules/convert.smk` (`guppy_tog`).

## Command line

```sh
phylobench-measure --out {prefix}.tog.run.json --threads {threads} \
    --tool-version 1.1.alpha23 --env-lock envs/pplacer.linux-64.pin.txt -- \
    guppy tog -o {prefix}.tog.tre {prefix}.jplace
```

stdout/stderr go to `{prefix}.tog.log`. On `status != ok` the rule touches the output tree so the job succeeds and `run.json` survives (`docs/harness.md` section 2). That DNF path is untested until T8 has a forced-DNF case.

## Version

- guppy ships in the `pplacer` conda package, 1.1.alpha23 (build `h01ede22_0`). `guppy --version` prints only `dev`; the proof is the package metadata (see `docs/audit.md` section 2).
- `--tool-version 1.1.alpha23` is hard-coded in the rule. If `envs/pplacer.yaml` is re-pinned, update the rule too (drift risk).

## Threads

guppy is single-threaded and has no thread flag. `threads` only sets the SLURM reservation and the harness `threads` field.

## Quirks

- Best placement only: `guppy tog` takes the top placement per query.
- A pquery with several names becomes a zero-length polytomy.
- The `{prefix}` file layout is provisional; T8 owns the layout.
- CI proof: `.github/workflows/convert.yml` runs `tests/test_guppy_tog.py` with the locked pplacer env (the test needs `GUPPY` and is skipped otherwise).

## License

GPL-3 (see `LICENSES.md`).
