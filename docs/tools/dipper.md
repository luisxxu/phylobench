# DIPPER

T15 adds the command line, threads, and quirks. This file currently covers only the output adapter (T25).

## Adapter: `phylobench.dipper_adapter.dipper_to_jplace`

`dipper_to_jplace(dipper_newick, backbone_newick, queries) -> dict` turns DIPPER's backbone-plus-queries Newick into a jplace doc (version 3) that `phylobench.accuracy.score_jplace` scores like any other arm.

- The doc is NaN-bearing: pendant and distal lengths are NaN. Score it in memory only; never write it for guppy or any strict-JSON consumer.
- Edges are matched by canonical leaf-hash splits, so DIPPER's rooting does not matter. It fails loudly (`ValueError`) when DIPPER's tree restricted to the backbone differs from the backbone (RF check), when a backbone leaf is missing, when an unrequested leaf is present, on duplicate queries, and when a query is named like a backbone leaf.
- Joint placement is documented, not fixed: DIPPER places queries together, so a query-only clade maps to the backbone edge where the clade attaches, and all its members get that edge.
- A query attached at a backbone vertex is placed on the edge above that vertex (the first child edge for the root vertex). Against a truth on a child edge of that vertex this scores edge distance 1 (pessimistic by up to 1).
- A query absent from DIPPER's output gets no record, so `score_jplace` returns a NaN row for it. Count NaN as missing; do not `skipna`.
- Cost: about 30 s for 100k queries on a 330k-leaf backbone (`tests/test_dipper_adapter_perf.py`, limit 60 s).
- UNVERIFIED: whether `dipper --add` preserves the backbone topology. T15 verifies it on fixture output; until then the RF check is the guard.
