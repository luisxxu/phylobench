# sepp-pinned (`align+place`, arm `sepp-pinned`)

SEPP in a pinned conda env, run on the shared Bfull\*/B10k\* backbones. Rule: `workflow/rules/sepp-pinned.smk` (`sepp_pinned_place`); arm file `config/arms/sepp-pinned.yaml`. Unaligned arm: it runs once per region with `aln=raw` on `results/{dataset}/regions/{region}/queries.fasta`.

## Command line

```sh
work={cell}/work; mkdir -p $work/out $work/tmp
phylobench-measure --out {cell}/align+place.run.json --threads {threads} \
    --tool-version 4.5.6 --env-lock envs/sepp.linux-64.pin.txt -- \
    run_sepp.py -t data/{dataset}/{backbone}/tree.nwk \
    -r data/{dataset}/{backbone}/raxml-info.txt \
    -a data/{dataset}/{backbone}/aligned-dna-sequences.fasta \
    -f results/{dataset}/regions/{region}/queries.fasta \
    -x {threads} -A 1000 -P 5000 -o sepp -d $work/out -p $work/tmp
python workflow/scripts/sepp_jplace.py $work/out/sepp_placement.json {cell}/placements.jplace; rm -rf $work
```

stdout/stderr go to `{cell}/align+place.log`. On `status != ok` the rule touches `placements.jplace` so `run.json` survives; an `ok` status with no placement file hard-fails the job. The `convert` phase (`guppy tog`) follows in the same cell via the shared `guppy_tog` rule (`docs/tools/guppy.md`).

## Version

`envs/sepp.yaml`, lock `envs/sepp.linux-64.pin.txt`: sepp 4.5.6 (`py313h106432d_5`), pplacer 1.1.alpha23 (`h01ede22_0`), hmmer 3.4, openjdk 25.0.2 (SEPP runs `java -jar seppJsonMerger.jar` and bioconda sepp does not declare java). `--tool-version 4.5.6` is hard-coded in the rule; update it with the pin.

## Threads

`-x {threads}` (config `threads`, 16 in production). Each SEPP worker runs `pplacer -j 1` and `hmmer --cpu 1`.

## License

SEPP is GPL-3. The engines executed are pplacer (GPL-3), HMMER (BSD-3) and the OpenJDK-run JSON merger; `guppy tog` runs in the `convert` phase (see `LICENSES.md`).

## Phases

`align+place` (SEPP aligns with HMMER and places with pplacer in one call, so it is one measured phase) and `convert`. T23 sums `align+place` and `convert` to compare against `sepp-prod`.

## Quirks

- `run_sepp.py` is used instead of `run-sepp.sh` (a deviation from the design section 2 wording) so the shared `guppy_tog` convert rule can consume the jplace.
- SEPP writes jplace v1 with `[n]` edge numbers; T26's `workflow/scripts/sepp_jplace.py` (`phylobench.jplace.from_sepp`) rewrites it to v3 after an `ok` run, outside the harness (unmeasured), as for `sepp-prod`.
- `resolve_polytomies`: SEPP resolves the trifurcating backbone root into an extra edge (both root-pair edges map to one backbone edge in scoring). A non-root polytomy hard-fails scoring.
- Internal labels get randomised renaming; `rename-json.py` is ignored because scoring is by leaf split and leaves are never renamed.
- Config comes from the dist-info `home.path` to `$PREFIX/share/sepp/sepp/main.config`, not `~/.sepp`; CI checks this in `.github/workflows/envs.yml` (the sepp check).
- Outdir prefix collision: SEPP output names depend on `-o` and the output dir, hence a fresh cell-local work dir per job.
- `-A`/`-P` are clamped to the number of taxa, so the smoke run is a single subset.
- Seed 297834.
- Output is named `{prefix}_placement.json`.
- Scratch is cell-local because a tmpfs `TMPDIR` would count toward cgroup `memory.peak`.
- `smoke_min_exact` is not set yet; add it once CI reports the exact count.

## Provenance

The pplacer `url#md5` line is identical across the pplacer, sepp and bscampp locks (the lock-guard in `.github/workflows/envs.yml` enforces it). Every `run.json` records `env_lock_sha256`.
