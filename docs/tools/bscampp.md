# bscampp family (`place`, arms `pplacer-scampp`, `epa-ng-scampp`, `bscampp`)

Three arms from one pinned `bscampp` install (`envs/bscampp.yaml`, lock `envs/bscampp.linux-64.pin.txt`, wheel and config written by `envs/bscampp.post-deploy.sh`). Rules: `workflow/rules/bscampp.smk` (`pplacer_scampp_place`, `epa_ng_scampp_place`, `bscampp_place`); arm files `config/arms/{pplacer-scampp,epa-ng-scampp,bscampp}.yaml`. Same cell layout as `docs/tools/pplacer.md`; the `convert` phase (`guppy tog`) follows (`docs/tools/guppy.md`).

## Command lines

`{cell}` is `results/{dataset}/{backbone}/{region}/{aln}/<arm>`. Each is preceded by `rm -rf {cell}/tmp0` and run as `OMP_NUM_THREADS={threads} phylobench-measure --out {cell}/place.run.json --threads {threads} --tool-version '<version>' --env-lock envs/bscampp.linux-64.pin.txt -- ...` with output to `{cell}/place.log`.

```sh
# bscampp            --tool-version 'bscampp 1.0.8 (epa-ng 0.3.8)'
run_bscampp.py --placement-method epa-ng -i results/ref-prep/{dataset}/{backbone}/epa-ng.model \
    -t {tree} -a {ref} -q {queries} -d {cell} -o placements.jplace \
    --threads {threads} --cpus-per-job 2 -b 2000 -V 5

# epa-ng-scampp      --tool-version 'bscampp 1.0.8 (epa-ng 0.3.8)'
run_scampp.py --placement-method epa-ng -i results/ref-prep/{dataset}/{backbone}/epa-ng.model \
    -t {tree} -a {ref} -q {queries} -d {cell} -o placements.jplace \
    --threads {threads} --cpus-per-job 1 -b 2000

# pplacer-scampp     --tool-version 'bscampp 1.0.8 (pplacer 1.1.alpha23)'
run_scampp.py --placement-method pplacer -i data/{dataset}/{backbone}/raxml-info.txt \
    -t {tree} -a {ref} -q {queries} -d {cell} -o placements.jplace \
    --threads {threads} --cpus-per-job 1 -b 2000
```

`{tree}`, `{ref}`, `{queries}` are the same inputs as the `pplacer` arm. Subtree size 2000 and, for `bscampp`, votes 5 are passed explicitly so they are recorded in the command.

## Versions

- bscampp 1.0.8, installed by wheel hash in `envs/bscampp.post-deploy.sh`. The wheel is not part of the conda lock, so `env_lock_sha256` does not cover it.
- pplacer 1.1.alpha23 and epa-ng 0.3.8, the same lock lines as `envs/pplacer.linux-64.pin.txt` and `envs/epa-ng.linux-64.pin.txt`.
- The `tool` column is the package name (`bscampp`); the arm name distinguishes the three.

## Threads

`--threads` is the worker budget. SCAMPP runs a pool of `num_cpus` workers, each given `cpus_per_job` threads, so SCAMPP uses `--cpus-per-job 1` to keep the total at `{threads}`. BSCAMPP runs `threads // cpus_per_job` workers, so `--cpus-per-job 2` gives 8 workers at 16 threads. `OMP_NUM_THREADS={threads}` is set for the harness and its children.

## License

The wrapper is MIT. Engines actually executed: `pplacer-scampp` runs pplacer (GPL-3) and taxit (GPL-3); `epa-ng-scampp` and `bscampp` run EPA-ng (AGPL-3). All three run the bundled homology helper (`tools/hamming_distance/*`, MIT). See `LICENSES.md`.

## Config mechanism

`main.config` is written at env deploy to `$CONDA_PREFIX/share/bscampp/main.config` and `home.path` points there, so no run touches `~/.bscampp`. Each rule checks `place.log` (the config dump bscampp prints is the provenance of the binary paths):

- on `status=ok`: `main.config: $CONDA_PREFIX/share/bscampp/main.config` and `Configs.epang_path: $CONDA_PREFIX/bin/epa-ng` (EPA-ng arms) or `Configs.pplacer_path` and `Configs.taxit_path` under `$CONDA_PREFIX/bin` (`pplacer-scampp`);
- on every status: the log must not contain `Cannot find home.path`, `Cannot find main.config` or `Initializing the config file`, which bscampp prints when it regenerates `~/.bscampp`.

## Quirks

- No `--tmpdir`: bscampp writes `tmp0` inside the cell, which the rule removes first so reruns start clean.
- `-i` is the model for the engine: the T6 `epa-ng.model` string for EPA-ng, the RAxML-8 info file for pplacer. The path must not contain `bestModel`, which switches bscampp's parsing.
- bscampp calls `resolve_polytomies` on the backbone and may mis-handle a non-root polytomy. GG2 backbones are at risk; check in T18/T22.
- Queries with no homology to the backbone are dropped by bscampp.
- The fixture backbone has 100 leaves, below the subtree size, so the subtree is the whole tree.
- An engine OOM inside a pool worker surfaces as `status=error`, not an OOM kill of the harness child.
- Aligned arm: places the `hmmalign` or `oracle` aligned queries.
