import time

import numpy as np
from tree_helpers import coalescent_tree

from phylobench.dipper_adapter import dipper_to_jplace
from phylobench.jplace import to_newick

TIME_LIMIT_S = 60


def test_330k_backbone_with_100k_queries():
    rng = np.random.default_rng(0)
    full = coalescent_tree(430_000, rng)
    labels = sorted(leaf.label for leaf in full.traverse_leaves())
    queries = [str(q) for q in rng.choice(labels, size=100_000, replace=False)]
    backbone = full.extract_tree_without(set(queries))
    dipper_newick, backbone_newick = to_newick(full), to_newick(backbone)

    start = time.perf_counter()
    doc = dipper_to_jplace(dipper_newick, backbone_newick, queries)
    elapsed = time.perf_counter() - start

    print(f"330k-backbone 100k-query dipper_to_jplace: {elapsed:.1f} s")
    assert elapsed < TIME_LIMIT_S
    assert len(doc["placements"]) == len(queries)
    assert all(rec["p"][0][0] is not None for rec in doc["placements"])
