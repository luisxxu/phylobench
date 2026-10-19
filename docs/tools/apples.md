# APPLES-2 (`place`, arm `apples2`)

Distance-based placement arm. Rule: `workflow/rules/apples.smk` (`apples2_place`); arm file `config/arms/apples2.yaml`. It covers every `alns` entry, including `hmmalign` and `oracle`.

## Command line

```sh
phylobench-measure --out {cell}/place.run.json --threads {threads} \
    --tool-version 2.0.11 --env-lock envs/apples.linux-64.pin.txt -- \
    run_apples.py -s data/{dataset}/{backbone}/aligned-dna-sequences.fasta \
    -q results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta \
    -t results/ref-prep/{dataset}/{backbone}/me-tree.nwk \
    -T {threads} -D -o {cell}/placements.jplace
```

stdout/stderr go to `{cell}/place.log`. On `status != ok` the rule touches `placements.jplace` so `run.json` survives (`docs/harness.md` section 2). The `convert` phase (`guppy tog`) follows in the same cell (`docs/tools/guppy.md`).

## Inputs

- `-s`: backbone alignment; `-q`: queries only (the `hmmalign` or `oracle` aligned queries).
- `-t`: the T6 `me-tree.nwk`, built once per backbone with FastTree `-nome -noml` (`docs/tools/fasttree.md`).
- `-D` is required. Without it the bundled `tools/FastTree-linux` re-estimates the tree per job.

## Version

`apples==2.0.11` is installed from PyPI by sha256 in `envs/apples.post-deploy.sh`, which Snakemake runs after the pin-file env is created (needs network at env creation). `--env-lock` hashes the conda pin file only. `--tool-version 2.0.11` is hard-coded in the rule; update it with the pin.

## Threads

`-T {threads}` is always passed (default 0 = all cores), plus the SLURM reservation and the harness `threads` field.

## Output and scoring

v3 jplace whose embedded tree is `me-tree.nwk` with APPLES postorder `{n}` numbering, not `tree.nwk`'s. `evaluate` maps edges by split from the embedded tree (`jplace_edge_nodes` in `src/phylobench/accuracy.py`), so rooting and numbering differences are fine. This is looser than the "edge numbers as in tree.nwk" wording in `docs/harness.md` section 7.

The `likelihood` field is really the least-squares error, and `like_weight_ratio` is always 1.

## License

GPL-3. Engines executed: APPLES and TreeCluster 1.0.5 (GPL-3). The bundled FastTree is not run (`-D`); ME branch lengths come from the pinned FastTree 2.2.0 in `ref-prep`. See `LICENSES.md`.

## Quirks

- Truth caveat: GG2 16S truth was seeded by uDance (APPLES-2 plus ML refinement), so apples2 results on GG2 truth are caveated and may be favoured. Do not read an apples2 win as a placement-quality win without the T16 simulated cross-check.
- Known issue: upstream `jutil.join_jplace` keeps the first query record even if it is unplaceable (edge -1), which makes `evaluate` raise `KeyError`. Not reachable on the fixture; flagged for T22.
