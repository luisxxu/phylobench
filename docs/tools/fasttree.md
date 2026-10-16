# FastTree 2.2.0

Env: `envs/fasttree.yaml` (lock `envs/fasttree.linux-64.pin.txt`).

## ME branch lengths for APPLES-2 (`ref-prep` phase)

Rule `refprep_me_tree` re-estimates minimum-evolution branch lengths on the fixed backbone topology:

```
OMP_NUM_THREADS={threads} FastTreeMP -nosupport -nt -nome -noml -intree {tree} {aln} > {out}
```

- Run through `phylobench-measure`; outputs are `results/ref-prep/{dataset}/{backbone}/me-tree.nwk` and `me-tree.run.json`, consumed by T8, T10, and T14. DNF handling follows `docs/harness.md` section 2.
- `FastTreeMP` takes its thread count from `OMP_NUM_THREADS` (set from the rule's threads).
- `-gtr` is omitted: with `-noml` the branch lengths are ME-only, and the fixed model is not used.
- Whether the output is identical across thread counts is unverified.
- Input: `data/{dataset}/{backbone}/{tree.nwk,aligned-dna-sequences.fasta}`. The backbone must be strictly bifurcating with no duplicate aligned rows; T16/T18 must guarantee this.
