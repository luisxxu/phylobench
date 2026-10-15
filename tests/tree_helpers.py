import treeswift

from phylobench.jplace import edge_numbered_tree


def to_newick(tree: treeswift.Tree, numbered: bool = False) -> str:
    """Iterative Newick writer; `numbered` adds a preorder `{n}` to non-root nodes."""
    out: list[str] = []
    stack: list[list] = []
    counter = 0

    def tail(node, number):
        length = "" if node.edge_length is None else f":{node.edge_length!r}"
        edge = "" if number is None else f"{{{number}}}"
        return f"{node.label or ''}{length}{edge}"

    def enter(node, is_root):
        nonlocal counter
        number = None
        if numbered and not is_root:
            number, counter = counter, counter + 1
        if node.children:
            out.append("(")
            stack.append([node, iter(node.children), number, True])
        else:
            out.append(tail(node, number))

    enter(tree.root, True)
    while stack:
        frame = stack[-1]
        child = next(frame[1], None)
        if child is None:
            out.append(")" + tail(frame[0], frame[2]))
            stack.pop()
            continue
        if not frame[3]:
            out.append(",")
        frame[3] = False
        enter(child, False)
    return "".join(out) + ";"


def numbered_doc(backbone_newick: str, fields=None) -> dict:
    """Empty jplace dict whose tree is the backbone numbered in preorder."""
    tree = treeswift.read_tree_newick(backbone_newick)
    fields = fields or [
        "edge_num",
        "like_weight_ratio",
        "distal_length",
        "pendant_length",
    ]
    return {
        "version": 3,
        "fields": fields,
        "tree": to_newick(tree, numbered=True),
        "placements": [],
    }


def edge_numbers(doc: dict) -> dict[str, int]:
    """Node label -> edge number for every labelled node of the doc's tree."""
    _, nodes = edge_numbered_tree(doc["tree"])
    return {node.label: number for number, node in nodes.items() if node.label}


def place(doc: dict, query: str, edge_num: int, lwr: float = 1.0) -> None:
    doc["placements"].append({"p": [[edge_num, lwr, 0.0, 0.1]], "n": [query]})
