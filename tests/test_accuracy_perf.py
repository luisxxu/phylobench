import time

import numpy as np
import treeswift
from tree_helpers import to_newick

from phylobench.accuracy import score_jplace

TIME_LIMIT_S = 60


def coalescent_tree(leaf_count: int, rng, zero_fraction: float = 0.1) -> treeswift.Tree:
    active = [treeswift.Node(label=f"t{i}") for i in range(leaf_count)]
    while len(active) > 1:
        picks = []
        for _ in range(2):
            i = int(rng.integers(len(active)))
            active[i], active[-1] = active[-1], active[i]
            picks.append(active.pop())
        parent = treeswift.Node()
        for child in picks:
            parent.add_child(child)
        active.append(parent)
    tree = treeswift.Tree()
    tree.root = active[0]
    for node in tree.traverse_preorder():
        if node is not tree.root:
            node.edge_length = (
                0.0 if rng.random() < zero_fraction else float(rng.exponential(0.1))
            )
    return tree


def time_scoring(
    backbone_newick: str, reference_newick: str, queries: list[str], rng
) -> float:
    numbered = to_newick(treeswift.read_tree_newick(backbone_newick), numbered=True)
    edges = rng.integers(numbered.count("{"), size=len(queries))
    doc = {
        "version": 3,
        "fields": ["edge_num", "like_weight_ratio", "distal_length", "pendant_length"],
        "tree": numbered,
        "placements": [
            {"p": [[int(e), 1.0, 0.0, 0.1]], "n": [q]}
            for q, e in zip(queries, edges, strict=True)
        ],
    }
    start = time.perf_counter()
    rows = score_jplace(doc, backbone_newick, reference_newick, queries)
    elapsed = time.perf_counter() - start
    assert len(rows) == len(queries)
    assert not np.isnan([r["edge_dist"] for r in rows]).any()
    return elapsed


def test_331k_leaf_coalescent_tree_with_zero_length_edges():
    rng = np.random.default_rng(0)
    reference = coalescent_tree(331_000, rng)
    labels = sorted(leaf.label for leaf in reference.traverse_leaves())
    queries = [str(q) for q in rng.choice(labels, size=1000, replace=False)]
    backbone = reference.extract_tree_without(set(queries))

    elapsed = time_scoring(to_newick(backbone), to_newick(reference), queries, rng)
    print(f"331k-leaf coalescent score_jplace: {elapsed:.1f} s")
    assert elapsed < TIME_LIMIT_S


def caterpillar(labels: list[str]) -> str:
    opening = "".join(f"({label}:1," for label in labels[:-1])
    return f"{opening}{labels[-1]}:1{'):1' * (len(labels) - 2)});"


def test_50k_leaf_caterpillar():
    rng = np.random.default_rng(1)
    labels = [f"t{i}" for i in range(51_000)]
    queries = [str(q) for q in rng.choice(labels[1:-1], size=1000, replace=False)]
    held_out = set(queries)

    backbone = caterpillar([x for x in labels if x not in held_out])
    elapsed = time_scoring(backbone, caterpillar(labels), queries, rng)
    print(f"50k-leaf caterpillar score_jplace: {elapsed:.1f} s")
    assert elapsed < TIME_LIMIT_S
