"""Minimal jplace reader and writer: top-LWR placement, edge-numbered trees."""

import json
import math
import re
from dataclasses import dataclass

import treeswift

REQUIRED_FIELDS = ("edge_num", "like_weight_ratio", "distal_length", "pendant_length")

# A node's edge number trails its branch length: `name:0.1{3}`. Moving it into the
# label (`name@3:0.1`) lets a plain Newick parser carry it.
_EDGE_NUMBER = re.compile(r"(:[^,()\[\]{};]*)?\{(\d+)\}")
_NUMBERED_LABEL = re.compile(r"(.*)@(\d+)")
# pplacer prints floats with C printf (`-nan`, `inf`), which strict JSON rejects.
_NON_FINITE = re.compile(r"(?<=[\[,\s])(-?)(nan|inf)(?=[\],\s])")


@dataclass(frozen=True)
class Placement:
    query: str
    edge_num: int
    like_weight_ratio: float
    distal_length: float
    pendant_length: float


def parse_jplace(text: str) -> dict:
    """json.loads that also reads pplacer's bare `nan` / `inf` numbers."""
    return json.loads(
        _NON_FINITE.sub(lambda m: "NaN" if m[2] == "nan" else f"{m[1]}Infinity", text)
    )


def top_placements(doc: dict) -> list[Placement]:
    """Highest-LWR placement per query; ties keep the first row (arbitrary but fixed).

    Rows with a non-finite LWR are ignored (pplacer's likelihoods underflow on very
    large backbones), and a record left with no rows is skipped, so its queries are
    unplaced.
    """
    if doc.get("version") != 3:
        raise ValueError(f"unsupported jplace version: {doc.get('version')!r}")
    fields = doc["fields"]
    missing = [f for f in REQUIRED_FIELDS if f not in fields]
    if missing:
        raise ValueError(f"jplace fields missing: {missing}")
    edge, lwr, distal, pendant = (fields.index(f) for f in REQUIRED_FIELDS)

    placements: dict[str, Placement] = {}
    for record in doc["placements"]:
        rows = [r for r in record["p"] if math.isfinite(r[lwr])]
        if not rows:
            continue
        row = max(rows, key=lambda r: r[lwr])
        names = record["n"] if "n" in record else [name for name, _ in record["nm"]]
        for name in names:
            if name in placements:
                raise ValueError(f"duplicate query in jplace: {name}")
            placements[name] = Placement(
                name,
                int(row[edge]),
                float(row[lwr]),
                float(row[distal]),
                float(row[pendant]),
            )
    return list(placements.values())


def edge_numbered_tree(newick: str) -> tuple[treeswift.Tree, dict[int, treeswift.Node]]:
    """Parse a jplace `{edge_num}` tree; unnamed nodes get label None."""
    tree = treeswift.read_tree_newick(_EDGE_NUMBER.sub(r"@\2\1", newick))
    nodes: dict[int, treeswift.Node] = {}
    for node in tree.traverse_preorder():
        match = node.label and _NUMBERED_LABEL.fullmatch(node.label)
        if not match:
            continue
        name, number = match.groups()
        node.label = name or None
        if int(number) in nodes:
            raise ValueError(f"duplicate edge number: {number}")
        nodes[int(number)] = node
    return tree, nodes


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


_SEPP_EDGE_NUMBER = re.compile(r"(:[^,()\[\]{};]*)\[(\d+)\]")


def from_sepp(doc: dict) -> dict:
    """SEPP's jplace v1 (`:len[N]` edges) as v3 (`:len{N}`); everything else is kept."""
    if doc.get("version") != 1:
        raise ValueError(f"expected jplace version 1, got {doc.get('version')!r}")
    missing = [f for f in REQUIRED_FIELDS if f not in doc["fields"]]
    if missing:
        raise ValueError(f"jplace fields missing: {missing}")
    if "{" in doc["tree"]:
        raise ValueError("jplace tree already has {edge} numbers")
    tree, count = _SEPP_EDGE_NUMBER.subn(r"\1{\2}", doc["tree"])
    if not count:
        raise ValueError("no [edge] numbers found in jplace tree")
    return {**doc, "version": 3, "tree": tree}
