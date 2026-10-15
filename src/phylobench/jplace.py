"""Minimal jplace reader: top-LWR placement per query and edge-numbered trees."""

import re
from dataclasses import dataclass

import treeswift

REQUIRED_FIELDS = ("edge_num", "like_weight_ratio", "distal_length", "pendant_length")

# A node's edge number trails its branch length: `name:0.1{3}`. Moving it into the
# label (`name@3:0.1`) lets a plain Newick parser carry it.
_EDGE_NUMBER = re.compile(r"(:[^,()\[\]{};]*)?\{(\d+)\}")
_NUMBERED_LABEL = re.compile(r"(.*)@(\d+)")


@dataclass(frozen=True)
class Placement:
    query: str
    edge_num: int
    like_weight_ratio: float
    distal_length: float
    pendant_length: float


def top_placements(doc: dict) -> list[Placement]:
    """Highest-LWR placement per query; ties keep the first row (arbitrary but fixed).

    Records with an empty placement list are skipped, so those queries are unplaced.
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
        if not record["p"]:
            continue
        row = max(record["p"], key=lambda r: r[lwr])
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
