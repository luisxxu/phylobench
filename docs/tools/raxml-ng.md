# raxml-ng 2.0.3

Env: `envs/raxml-ng.yaml` (lock `envs/raxml-ng.linux-64.pin.txt`; the bioconda build pulls openmpi).

## Role

CI parse check only: `envs.yml` runs `raxml-ng --check` on the generated model files to prove the model string is valid. raxml-ng is not in any arm's run path.

## EPA-ng model file

EPA-ng cannot read the RAxML-8 GTRCAT `raxml-info.txt` shipped with the GG2 reference, so T6 converts it (`src/phylobench/models.py`, `workflow/scripts/write_epa_ng_model.py`). Rule `refprep_epa_ng_model` writes `results/ref-prep/{dataset}/{backbone}/epa-ng.model`, one line:

```
GTR{ac/ag/at/cg/ct/gt}+FU{A/C/G/T}+G4m{alpha}, noname = 1-{alignment width}
```

- Rates are read positionally and the info file's labels must be exactly `ac ag at cg ct gt`; anything else is an error.
- Parsing is strict: one partition, 6 rates, 4 frequencies, all finite and positive.
- `alpha` is passed through as shipped (1.0 for GTRCAT info files); models are never re-estimated.
- `FU` uses the shipped base frequencies as fixed values. pplacer instead defaults to empirical frequencies of the reference alignment (`--model-freqs` disables that), so the two arms differ slightly here.
- `G4m` is the mean-of-category discretisation, like pplacer.
- The width comes from the first record of the backbone alignment.

## Fixtures

`tests/fixtures/gg2-2022.10/raxml-info.txt` is the byte-identical `raxml-info.txt` member of the GG2 2022.10 SEPP reference (provenance and checksums in `tests/test_models.py`). The tiny fixture's info file is `tests/fixtures/tiny/raxml-info.txt`.

## Input contract (shared with T18)

`data/{dataset}/{backbone}/{tree.nwk,aligned-dna-sequences.fasta,raxml-info.txt}`.
