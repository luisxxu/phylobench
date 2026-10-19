# EPA-ng (`place`, arm `epa-ng`)

Control arm: plain EPA-ng without SCAMPP subtree selection. Rule: `workflow/rules/epa-ng.smk` (`epa_ng_place`); arm file `config/arms/epa-ng.yaml`.

## Command line

```sh
phylobench-measure --out {cell}/place.run.json --threads {threads} \
    --tool-version 'epa-ng 0.3.8' --env-lock envs/epa-ng.linux-64.pin.txt -- \
    epa-ng --ref-msa {ref} --tree {tree} --query {queries} \
    --model results/ref-prep/{dataset}/{backbone}/epa-ng.model \
    --outdir {cell} --threads {threads} --redo
```

Output goes to `{cell}/place.log`. EPA-ng always writes `epa_result.jplace`, so on `status=ok` the rule renames it to `placements.jplace`; on `status != ok` it touches `placements.jplace` instead, as the `pplacer` rule does.

## Version

`epa-ng=0.3.8` (`envs/epa-ng.yaml`, lock `envs/epa-ng.linux-64.pin.txt`).

## Threads

`--threads {threads}` (16 in production) plus the SLURM reservation and the harness `threads` field.

## License

AGPL-3, and the engine executed is EPA-ng itself (see `LICENSES.md`).

## Quirks

- `--model` is the T6 model string (`GTR{..}+FU{..}+G4m{1.0}`) with the shipped numbers, never re-estimated.
- B10k\* only (`design.md` section 5). T12 does not restrict it: arm config cannot name backbones, so T22 drops it at Bfull\* (see `TICKETS.md` T22).
- `--redo` lets a rerun overwrite the previous output in the cell.
