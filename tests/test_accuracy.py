import math

import numpy as np
import pytest
from tree_helpers import edge_numbers, place

from phylobench.accuracy import (
    Backbone,
    canonical,
    edge_distance,
    jplace_edge_nodes,
    leaf_hashes,
    placed_nodes,
    prepare_backbone,
    score_jplace,
    true_positions,
)
from phylobench.jplace import numbered_doc

# Binary except for the trifurcating root R; internal nodes are named for placement.
BINARY = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"


def run(backbone, reference, queries, placed, jplace_tree=None):
    doc = numbered_doc(jplace_tree or backbone)
    numbers = edge_numbers(doc)
    for query, label in placed.items():
        place(doc, query, numbers[label])
    return {
        row["query"]: row for row in score_jplace(doc, backbone, reference, queries)
    }


def node_of(backbone: Backbone, *leaves: str) -> int:
    sig = 0
    for leaf in leaves:
        sig ^= backbone.hashes[leaf]
    return backbone.key_to_node[canonical(sig, backbone.total)]


def test_leaf_hashes_do_not_depend_on_label_order():
    assert leaf_hashes(["a", "b", "c"]) == leaf_hashes(["c", "a", "b"])
    assert leaf_hashes(["a", "b"], seed=1) != leaf_hashes(["a", "b"], seed=2)


@pytest.mark.parametrize(
    ("placed_at", "expected"),
    [("a", 0), ("b", 1), ("X", 1), ("Z", 2), ("c", 4), ("f", 2)],
)
def test_binary_distances_from_mid_edge_truth(placed_at, expected):
    reference = "(((q:1,a:0.5)m:0.5,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    row = run(BINARY, reference, ["q"], {"q": placed_at})["q"]
    assert row["edge_dist"] == expected
    assert row["exact"] == float(expected == 0)


def test_delta_error_counts_missing_bipartitions():
    reference = "(((q:1,a:0.5)m:0.5,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    assert run(BINARY, reference, ["q"], {"q": "a"})["q"]["delta_err"] == 0
    assert run(BINARY, reference, ["q"], {"q": "b"})["q"]["delta_err"] == 1


def test_vertex_truth_is_zero_from_every_incident_edge():
    reference = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1,q:1)R;"
    rows = [
        run(BINARY, reference, ["q"], {"q": label})["q"]
        for label in ["X", "Z", "f", "a", "Y"]
    ]
    assert [r["edge_dist"] for r in rows] == [0, 0, 0, 1, 1]
    assert [r["delta_err"] for r in rows] == [0, 0, 0, 1, 1]


def test_internal_vertex_truth():
    reference = "((a:1,b:1,q:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    assert run(BINARY, reference, ["q"], {"q": "a"})["q"]["edge_dist"] == 0
    assert run(BINARY, reference, ["q"], {"q": "X"})["q"]["edge_dist"] == 0
    assert run(BINARY, reference, ["q"], {"q": "Z"})["q"]["edge_dist"] == 1


@pytest.mark.parametrize(
    "reference",
    [
        "(a:1,(b:1,q:1,(((c:1,d:1)Y:1,e:1)Z:1,f:1)R:1)X:1);",
        "(a:1,((((c:1,d:1)Y:1,e:1)Z:1,f:1)R:1,b:1,q:1)X:1);",
    ],
    ids=["carrier_below", "carrier_is_complement"],
)
def test_vertex_truth_when_reference_is_rooted_elsewhere(reference):
    dist = {
        label: run(BINARY, reference, ["q"], {"q": label})["q"]["edge_dist"]
        for label in ["a", "X", "f", "Y"]
    }
    assert dist == {"a": 0, "X": 0, "f": 1, "Y": 2}


def test_bifurcating_backbone_root_merges_edges_whichever_shape():
    internal_both = "((a:1,b:1)X:1,(c:1,d:1)Y:1)R;"
    reference = "((a:1,b:1)X:1,(c:1,d:1)Y:1,q:1)R;"
    for placed_at in ("X", "Y"):
        assert (
            run(internal_both, reference, ["q"], {"q": placed_at})["q"]["edge_dist"]
            == 0
        )
    assert run(internal_both, reference, ["q"], {"q": "a"})["q"]["edge_dist"] == 1

    leaf_and_internal = "((a:1,b:1)X:1,c:1)R;"
    reference = "((a:1,b:1)X:1,(q:1,c:0.5)m:0.5)R;"
    for placed_at in ("X", "c"):
        assert (
            run(leaf_and_internal, reference, ["q"], {"q": placed_at})["q"]["edge_dist"]
            == 0
        )


def test_merged_root_edge_weight():
    def merged_weight(newick):
        backbone = prepare_backbone(newick)
        return backbone.weight[node_of(backbone, "c", "d")]

    assert merged_weight("((a:1,b:1)X:0.0,(c:1,d:1)Y:0.5)R;") == 1
    assert merged_weight("((a:1,b:1)X:0.0,(c:1,d:1)Y:0.0)R;") == 0
    assert merged_weight("((a:1,b:1)X,(c:1,d:1)Y:0.0)R;") == 1


def test_none_length_edge_counts_as_nonzero():
    backbone = prepare_backbone("((a:1,b)X:0,c:1,d:1)R;")
    assert backbone.weight[node_of(backbone, "b")] == 1
    assert backbone.weight[node_of(backbone, "a", "b")] == 0


def test_separate_jplace_root_edges_map_to_one_node():
    backbone = prepare_backbone("((a:1,b:1)X:1,(c:1,d:1)Y:1)R;")
    doc = numbered_doc("((a:1,b:1)X:1,(c:1,d:1)Y:1)R;")
    mapping = jplace_edge_nodes(doc["tree"], backbone)
    numbers = edge_numbers(doc)
    assert mapping[numbers["X"]] == mapping[numbers["Y"]]


def test_differently_rooted_jplace_tree_maps():
    rerooted = "(a:1,(b:1,((c:1,d:1)Y:1,(e:1,f:1)Z:1)W:1)X:1)R;"
    backbone = "((a:1,b:1)X:1,(c:1,d:1)Y:1,(e:1,f:1)Z:1)R;"
    reference = "((a:1,b:1)X:1,(c:1,d:1)Y:1,(e:1,f:1)Z:1,q:1)R;"
    rows = run(backbone, reference, ["q"], {"q": "Y"}, jplace_tree=rerooted)
    assert rows["q"]["edge_dist"] == 0
    assert (
        run(backbone, reference, ["q"], {"q": "W"}, jplace_tree=rerooted)["q"][
            "edge_dist"
        ]
        == 0
    )


def test_duplicate_split_in_jplace_tree_raises():
    backbone = prepare_backbone(BINARY)
    unary = (
        "(((a:1{0},b:1{1})X:1{2})U:1{3},((c:1{4},d:1{5})Y:1{6},e:1{7})Z:1{8},f:1{9});"
    )
    with pytest.raises(ValueError, match="duplicates"):
        jplace_edge_nodes(unary, backbone)


def test_zero_length_internal_polytomy():
    backbone = "((a:1,b:1)X:0,(c:1,d:1)Y:0,e:1)R;"
    reference = "((a:1,b:1)X:0,(c:1,d:1)Y:0,e:1,q:1)R;"
    for label in ["a", "X", "Y", "e"]:
        assert run(backbone, reference, ["q"], {"q": label})["q"]["edge_dist"] == 0
    edge_truth = "(((q:1,a:0.5)m:0.5,b:1)X:0,(c:1,d:1)Y:0,e:1)R;"
    assert run(backbone, edge_truth, ["q"], {"q": "b"})["q"]["edge_dist"] == 1
    assert run(backbone, edge_truth, ["q"], {"q": "X"})["q"]["edge_dist"] == 0


def test_zero_length_true_edge_is_a_vertex():
    backbone = "((a:1,b:1)X:0,c:1,d:1)R;"
    reference = "(((a:1,b:1)X:0,q:1)m:0,c:1,d:1)R;"
    for label in ["a", "b", "c", "d", "X"]:
        assert run(backbone, reference, ["q"], {"q": label})["q"]["edge_dist"] == 0


def test_zero_length_leaf_edge():
    backbone = "((a:0,b:1)X:1,c:1,d:1)R;"
    reference = "(((q:1,a:0)m:0,b:1)X:1,c:1,d:1)R;"
    dist = {
        label: run(backbone, reference, ["q"], {"q": label})["q"]["edge_dist"]
        for label in "abXcd"
    }
    assert dist == {"a": 0, "b": 0, "X": 0, "c": 1, "d": 1}


def test_distance_across_zero_length_edges_is_collapsed():
    backbone = "((a:1,b:1)X:0,(c:1,d:1)Y:1,e:1)R;"
    reference = "(((q:1,a:0.5)m:0.5,b:1)X:0,(c:1,d:1)Y:1,e:1)R;"
    assert run(backbone, reference, ["q"], {"q": "Y"})["q"]["edge_dist"] == 1
    assert run(backbone, reference, ["q"], {"q": "c"})["q"]["edge_dist"] == 2


def test_mid_edge_attachment_does_not_break_split_check():
    reference = "((((q:1,a:0.5)m:0.5)n:0,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    assert run(BINARY, reference, ["q"], {"q": "a"})["q"]["exact"] == 1.0


def test_degree_two_nodes_in_reference_are_ignored():
    reference = "(((a:0.5)u:0.5,b:1)X:1,((c:1,d:1)Y:1,(q:1,e:0.5)m:0.5)Z:1,f:1)R;"
    assert run(BINARY, reference, ["q"], {"q": "e"})["q"]["exact"] == 1.0


def test_unary_backbone_node_raises():
    with pytest.raises(ValueError, match="duplicate"):
        prepare_backbone("(((a:1,b:1)X:1)U:1,c:1,d:1)R;")


def test_query_above_backbone_root_with_vertex_truth():
    reference = "(q:1,((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R:1)T;"
    assert run(BINARY, reference, ["q"], {"q": "X"})["q"]["edge_dist"] == 0
    assert run(BINARY, reference, ["q"], {"q": "a"})["q"]["edge_dist"] == 1


def test_query_above_backbone_root_with_edge_truth():
    backbone = "((a:1,b:1)X:1,(c:1,d:1)Y:1)R;"
    reference = "(q:1,(r:1,((a:1,b:1)X:1,(c:1,d:1)Y:1)R:1)S:1)T;"
    rows = run(backbone, reference, ["q", "r"], {"q": "X", "r": "Y"})
    assert rows["q"]["edge_dist"] == rows["r"]["edge_dist"] == 0


def test_clade_of_queries_shares_the_stem_edge():
    reference = (
        "(((a:1,((q:1,r:1)Q:1,s:1)P:1)m:0.5,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    )
    rows = run(BINARY, reference, ["q", "r", "s"], {"q": "a", "r": "b", "s": "a"})
    assert [rows[q]["edge_dist"] for q in "qrs"] == [0, 1, 0]


def test_singleton_cherry_of_queries_at_a_vertex():
    reference = "((a:1,b:1,(q:1,r:1)Q:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1)R;"
    rows = run(BINARY, reference, ["q", "r"], {"q": "X", "r": "a"})
    assert [rows[q]["edge_dist"] for q in "qr"] == [0, 0]


def test_reference_with_different_backbone_splits_raises():
    reference = "((a:1,c:1)X:1,((b:1,d:1)Y:1,e:1)Z:1,f:1,q:1)R;"
    with pytest.raises(ValueError, match="differs"):
        run(BINARY, reference, ["q"], {"q": "a"})


def test_reference_must_contain_backbone_and_queries():
    with pytest.raises(ValueError, match="backbone leaves and the queries"):
        true_positions(BINARY, ["q"], prepare_backbone(BINARY))


def test_extra_reference_leaves_leave_positions_unchanged():
    exact = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1,q:1)R;"
    extra = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1,q:1,spare:1)R;"
    backbone = prepare_backbone(BINARY)
    nodes, weights = true_positions(exact, ["q"], backbone)
    extra_nodes, extra_weights = true_positions(extra, ["q"], backbone)
    assert list(extra_nodes) == list(nodes)
    assert list(extra_weights) == list(weights)


def test_query_that_is_a_backbone_leaf_raises():
    with pytest.raises(ValueError, match="backbone leaves and the queries"):
        true_positions(BINARY, ["a"], prepare_backbone(BINARY))


def test_jplace_leaf_mismatch_raises():
    backbone = prepare_backbone("((a:1,b:1)X:1,c:1,d:1)R;")
    with pytest.raises(ValueError, match="leaves"):
        jplace_edge_nodes("((a:1{0},b:1{1}):1{2},c:1{3},z:1{4});", backbone)


def test_placed_name_not_in_queries_raises():
    reference = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1,q:1)R;"
    with pytest.raises(ValueError, match="not requested"):
        run(BINARY, reference, ["r"], {"q": "a"})


def test_duplicate_requested_query_raises():
    reference = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1,q:1)R;"
    with pytest.raises(ValueError, match="duplicate"):
        run(BINARY, reference, ["q", "q"], {"q": "X"})


def test_one_row_per_query_with_nan_for_unplaced():
    reference = "((a:1,b:1)X:1,((c:1,d:1)Y:1,e:1)Z:1,f:1,q:1,r:1)R;"
    rows = run(BINARY, reference, ["q", "r"], {"q": "X"})
    assert list(rows) == ["q", "r"]
    assert rows["q"]["edge_dist"] == 0.0
    assert all(math.isnan(rows["r"][k]) for k in ("edge_dist", "exact", "delta_err"))


def test_edge_distance_with_sepp_style_truth():
    backbone = prepare_backbone(BINARY)
    truth = np.array([node_of(backbone, "a")] * 3)
    placed = np.array(
        [node_of(backbone, "a"), node_of(backbone, "b"), node_of(backbone, "c")]
    )
    dist, exact, _ = edge_distance(backbone, placed, truth, backbone.weight[truth])
    assert dist.tolist() == [0, 1, 4]
    assert exact.tolist() == [True, False, False]


class RandomBackbone:
    """Small random backbone with a polytomous root and ~30% zero-length edges."""

    def __init__(self, rng):
        self.rng = rng
        self.children: dict[int, list[int]] = {}
        self.parent: dict[int, int] = {}
        self.length: dict[int, float] = {}
        self.name: dict[int, str] = {}
        leaves = [f"l{i}" for i in range(int(rng.integers(5, 12)))]
        self.root = self._grow(leaves, min_groups=3)

    def _new(self, name):
        node = len(self.name)
        self.name[node] = name
        self.children[node] = []
        return node

    def _grow(self, leaves, min_groups=2):
        node = self._new("")
        groups = int(self.rng.integers(min_groups, min(4, len(leaves)) + 1))
        cuts = sorted(
            self.rng.choice(np.arange(1, len(leaves)), size=groups - 1, replace=False)
        )
        for part in np.split(np.array(self.rng.permutation(leaves)), cuts):
            child = self._new(part[0]) if len(part) == 1 else self._grow(list(part))
            if len(part) > 1:
                self.name[child] = f"n{child}"
            self.children[node].append(child)
            self.parent[child] = node
            self.length[child] = (
                0.0 if self.rng.random() < 0.3 else float(self.rng.uniform(0.1, 1.0))
            )
        return node

    def edges(self):
        return list(self.parent)

    def internals(self):
        return [v for v in self.children if self.children[v]]

    def leaves_below(self, v):
        if not self.children[v]:
            return [self.name[v]]
        return [leaf for c in self.children[v] for leaf in self.leaves_below(c)]

    def newick(self, attachments=()):
        """attachments: (query, 'edge' | 'vertex', target, optional sister query)."""
        on_edge: dict[int, list] = {}
        at_vertex: dict[int, list] = {}
        for query, kind, target, extra in attachments:
            (on_edge if kind == "edge" else at_vertex).setdefault(target, []).append(
                (query, extra)
            )

        def pendant(query, extra):
            return f"({query}:1,{extra}:1)" if extra else f"{query}:1"

        def render(v):
            parts = [render(c) for c in self.children[v]] + [
                pendant(*a) for a in at_vertex.get(v, [])
            ]
            text = f"({','.join(parts)})" if parts else ""
            if v == self.root:
                return f"{text}R;"
            length = self.length[v]
            halves = length / 2
            attached = on_edge.get(v, [])
            if not attached:
                return f"{text}{self.name[v]}:{length}"
            text = f"{text}{self.name[v]}:{halves}"
            for n, (query, extra) in enumerate(attached):
                text = f"({text},{pendant(query, extra)})m{v}_{n}:{halves}"
            return text

        return render(self.root)


def bipartitions(tree):
    leaves = {n: frozenset([n.label]) for n in tree.traverse_leaves()}
    below = dict(leaves)
    for node in tree.traverse_postorder(leaves=False):
        below[node] = frozenset().union(*(below[c] for c in node.children))
    everything = below[tree.root]
    anchor = min(everything)
    length_of: dict[frozenset, float] = {}
    for node in tree.traverse_preorder():
        side = below[node]
        if node is tree.root or len(side) < 2 or len(side) > len(everything) - 2:
            continue
        side = everything - side if anchor in side else side
        length_of[side] = length_of.get(side, 0.0) + node.edge_length
    return {side for side, total in length_of.items() if total > 0}


class Oracle:
    """Collapsed-tree distances by explicit graph search."""

    def __init__(self, tree: RandomBackbone):
        self.tree = tree
        group = {v: v for v in tree.children}

        def find(v):
            while group[v] != v:
                v = group[v]
            return v

        for v, p in tree.parent.items():
            if tree.length[v] == 0:
                group[find(v)] = find(p)
        self.group = {v: find(v) for v in tree.children}
        neighbours: dict[int, set] = {g: set() for g in set(self.group.values())}
        for v, p in tree.parent.items():
            if tree.length[v] > 0:
                neighbours[self.group[v]].add(self.group[p])
                neighbours[self.group[p]].add(self.group[v])
        self.hops = {g: self._bfs(neighbours, g) for g in neighbours}

    @staticmethod
    def _bfs(neighbours, start):
        seen, frontier = {start: 0}, [start]
        while frontier:
            frontier = [
                n
                for g in frontier
                for n in neighbours[g]
                if n not in seen and not seen.update({n: seen[g] + 1})
            ]
        return seen

    def edge_item(self, v):
        if self.tree.length[v] == 0:
            return ("V", self.group[v])
        return ("E", frozenset([self.group[v], self.group[self.tree.parent[v]]]))

    def distance(self, x, y):
        if x[0] == "E" and y[0] == "E":
            if x == y:
                return 0
            return 1 + min(self.hops[a][b] for a in x[1] for b in y[1])
        if x[0] == "V" and y[0] == "V":
            return self.hops[x[1]][y[1]]
        edge, vertex = (x, y) if x[0] == "E" else (y, x)
        return min(self.hops[a][vertex[1]] for a in edge[1])


def random_attachments(tree: RandomBackbone, rng):
    attachments, truths = [], {}
    for i in range(int(rng.integers(1, 4))):
        query = f"q{i}"
        extra = f"x{i}" if rng.random() < 0.3 else None
        if rng.random() < 0.5:
            target = int(rng.choice(tree.edges()))
            attachments.append((query, "edge", target, extra))
            truths[query] = ("edge", target)
        else:
            target = int(rng.choice(tree.internals()))
            attachments.append((query, "vertex", target, extra))
            truths[query] = ("vertex", target)
    return attachments, truths


@pytest.mark.filterwarnings("ignore:Deleting a node:UserWarning")
@pytest.mark.parametrize("seed", range(200))
def test_edge_distance_and_delta_error_match_oracles(seed):
    import treeswift

    rng = np.random.default_rng(seed)
    tree = RandomBackbone(rng)
    backbone = prepare_backbone(tree.newick())
    oracle = Oracle(tree)

    def index_of(v):
        return 0 if v == tree.root else node_of(backbone, *tree.leaves_below(v))

    attachments, truths = random_attachments(tree, rng)
    queries = [a[0] for a in attachments] + [a[3] for a in attachments if a[3]]
    reference = tree.newick(attachments)
    true_node, true_weight = true_positions(reference, queries, backbone)

    for i, (query, kind, target, _) in enumerate(attachments):
        expected_index = index_of(target)
        expected_weight = backbone.weight[expected_index] if kind == "edge" else 0
        assert (true_node[i], true_weight[i]) == (expected_index, expected_weight)

        truth_item = (
            oracle.edge_item(target) if kind == "edge" else ("V", oracle.group[target])
        )
        placed = tree.edges()
        dist, exact, delta = edge_distance(
            backbone,
            np.array([index_of(p) for p in placed]),
            np.full(len(placed), true_node[i]),
            np.full(len(placed), true_weight[i]),
        )
        assert dist.tolist() == [
            oracle.distance(oracle.edge_item(p), truth_item) for p in placed
        ]
        assert exact.tolist() == [d == 0 for d in dist]

        restricted = treeswift.read_tree_newick(reference).extract_tree_without(
            set(queries) - {query}
        )
        truth_bips = bipartitions(restricted)
        for j in rng.choice(len(placed), size=min(5, len(placed)), replace=False):
            estimate = bipartitions(
                treeswift.read_tree_newick(
                    tree.newick([(query, "edge", placed[j], None)])
                )
            )
            assert delta[j] == len(truth_bips - estimate)


def test_placed_nodes_returns_rows_and_backbone_nodes_of_placed_queries():
    backbone = prepare_backbone(BINARY)
    doc = numbered_doc(BINARY)
    numbers = edge_numbers(doc)
    place(doc, "r", numbers["a"])
    place(doc, "p", numbers["c"])
    rows, nodes = placed_nodes(doc, backbone, ["p", "q", "r"])
    assert rows.tolist() == [2, 0]
    assert nodes.tolist() == [node_of(backbone, "a"), node_of(backbone, "c")]


def test_placed_nodes_unplaced_queries_are_absent():
    rows, nodes = placed_nodes(numbered_doc(BINARY), prepare_backbone(BINARY), ["q"])
    assert rows.size == 0 and nodes.size == 0


def test_placed_nodes_rejects_unknown_and_duplicate_queries():
    backbone = prepare_backbone(BINARY)
    doc = numbered_doc(BINARY)
    place(doc, "q", edge_numbers(doc)["a"])
    with pytest.raises(ValueError, match="not requested"):
        placed_nodes(doc, backbone, ["r"])
    with pytest.raises(ValueError, match="duplicate"):
        placed_nodes(doc, backbone, ["q", "q"])
