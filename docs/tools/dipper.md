# DIPPER

Version 0.1.5 is pinned (`envs/dipper.yaml`, `envs/dipper.linux-64.pin.txt`); the binary's `--version` prints 0.1.4. License: MIT.

## Arms and commands

Rules: `workflow/rules/dipper.smk`. Each measured `place` phase writes `dipper.nwk` and `place.run.json`; the unmeasured `dipper_to_jplace` rule then writes `placements.jplace` (skipped to an empty file when the gate `run.json` is not ok). `-k`/`-s` are accepted but ignored in 0.1.5 (fixed k=15, sketch 1000).

| Arm | Command | Input |
|---|---|---|
| `dipper-m` | `dipper_cpu --add -i m -d 2 -T {threads} -t tree.nwk -I in.fasta -O dipper.nwk` | backbone alignment + masked query alignment, concatenated (`awk 1`) |
| `dipper-m-gpu` | `apptainer exec --nv dipper-gpu.sif dipper --add -i m -d 2 -t tree.nwk -I in.fasta -O dipper.nwk` | same file as `dipper-m` |
| `dipper-r` | `dipper_cpu --add -i r -T {threads} -t tree.nwk -I in.fasta -O dipper.nwk` | backbone + raw queries, ungapped (`workflow/scripts/ungap_concat.awk`) |

## Threads and versions

- CPU: native `-T {threads}`. GPU: the build has no thread flag; the harness records 16 and SLURM bounds it.
- GPU version is the image digest in `envs/dipper-gpu.image.txt` (UNVERIFIED until T27). The image is built from upstream main, not the 0.1.5 tag. The `.sif` is pulled once to `data/images/dipper-gpu.sif` by the local rule `dipper_gpu_sif`.
- The GPU arm needs `apptainer` and an NVIDIA driver on the host (`--nv`). It is manual tier: cluster only, after T27, excluded from the CI smoke matrix.

## Quirks

- `-d` default mismatch: the CPU binary defaults to `-d 1` despite its help text, so `-d 2` is passed explicitly.
- The backbone must be single-line Newick.
- DIPPER randomly resolves polytomies, so backbones must be binary except for the root (pre-flight for T6/T21).
- `--no-shuffle` is a no-op in `--add`.
- The combined FASTA is a `temp()` file (about 330k rows per cell at Bfull\*), so it needs local temp disk.
- `dipper-r` runs at B10k\* only, selected by the run config (T21/T22), not by the arm file.
- The adapter step is unmeasured (it is evaluation plumbing, not placement cost).

## Adapter: `phylobench.dipper_adapter.dipper_to_jplace`

`dipper_to_jplace(dipper_newick, backbone_newick, queries) -> dict` turns DIPPER's backbone-plus-queries Newick into a jplace doc (version 3) that `phylobench.accuracy.score_jplace` scores like any other arm.

- The doc is NaN-bearing: pendant and distal lengths are NaN. `workflow/scripts/dipper_to_jplace.py` writes it to disk (bare `NaN`) for `evaluate` only. Dipper arms must never list a `convert` phase: `guppy_tog` matches any `*.jplace` and would choke on NaN.
- Edges are matched by canonical leaf-hash splits, so DIPPER's rooting does not matter. It fails loudly (`ValueError`) when DIPPER's tree restricted to the backbone differs from the backbone (RF check), when a backbone leaf is missing, when an unrequested leaf is present, on duplicate queries, and when a query is named like a backbone leaf.
- Joint placement is documented, not fixed: DIPPER places queries together, so a query-only clade maps to the backbone edge where the clade attaches, and all its members get that edge.
- A query attached at a backbone vertex is placed on the edge above that vertex (the first child edge for the root vertex). Against a truth on a child edge of that vertex this scores edge distance 1 (pessimistic by up to 1).
- A query absent from DIPPER's output gets no record, so `score_jplace` returns a NaN row for it. Count NaN as missing; do not `skipna`.
- Cost: about 30 s for 100k queries on a 330k-leaf backbone (`tests/test_dipper_adapter_perf.py`, limit 60 s).
- Topology preservation under `--add` is verified by the CI smoke (pending the first green run); the RF check is the guard.
