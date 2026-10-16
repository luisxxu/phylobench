import treeswift

from phylobench.jplace import edge_numbered_tree


def edge_numbers(doc: dict) -> dict[str, int]:
    """Node label -> edge number for every labelled node of the doc's tree."""
    _, nodes = edge_numbered_tree(doc["tree"])
    return {node.label: number for number, node in nodes.items() if node.label}


def place(doc: dict, query: str, edge_num: int, lwr: float = 1.0) -> None:
    doc["placements"].append({"p": [[edge_num, lwr, 0.0, 0.1]], "n": [query]})


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
