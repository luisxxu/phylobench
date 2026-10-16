"""Held-out clade/singleton selection and pruning of the tree and alignment.

Pure functions; the same stem-parent rule as simulate.choose_heldout keeps every group
on its own backbone bipartition. Taxonomy rank 0 is the highest rank in each lineage.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import treeswift

from phylobench.accuracy import contract_bifurcating_root

_MIXED = object()


@dataclass(frozen=True)
class HeldOutGroup:
    group_id: str
    group_type: str
    rank: int | None
    taxon: str
    members: tuple[str, ...]


def monophyletic_taxa(
    tree: treeswift.Tree, lineages: dict[str, tuple[str, ...]], rank: int
) -> dict[str, treeswift.Node]:
    """Taxa at `rank` whose leaves are exactly one node's leaves; value is that node.

    A leaf with no lineage entry at `rank` belongs to no taxon, so it breaks the
    monophyly of any node above it.
    """
    taxon_of = {}
    totals: dict[str, int] = {}
    for label, lineage in lineages.items():
        if len(lineage) > rank:
            taxon_of[label] = lineage[rank]
            totals[lineage[rank]] = totals.get(lineage[rank], 0) + 1

    common: dict[treeswift.Node, object] = {}
    size: dict[treeswift.Node, int] = {}
    for node in tree.traverse_postorder():
        if node.is_leaf():
            common[node] = taxon_of.get(node.label, _MIXED)
            size[node] = 1
            continue
        kinds = {common[c] for c in node.children}
        common[node] = kinds.pop() if len(kinds) == 1 else _MIXED
        size[node] = sum(size[c] for c in node.children)

    found: dict[str, treeswift.Node] = {}
    for node in tree.traverse_preorder():
        taxon = common[node]
        if taxon is not _MIXED and taxon not in found and size[node] == totals[taxon]:
            found[taxon] = node
    return found


def _is_free(node, has_heldout, stem_parents):
    parent = node.parent
    if parent is None or parent.parent is None or has_heldout[parent]:
        return False
    while parent is not None:
        if parent in stem_parents:
            return False
        parent = parent.parent
    return True


def choose_groups(
    tree: treeswift.Tree,
    lineages: dict[str, tuple[str, ...]],
    rank_counts: dict[int, int],
    n_singletons: int,
    max_leaves: int,
    rng: np.random.Generator,
    min_leaves: int = 25,
) -> list[HeldOutGroup]:
    """Disjoint monophyletic clades (highest rank first) then singletons.

    A group is eligible only if its stem parent is not the root and the parent's
    subtree holds no other held-out leaf, and no earlier stem parent lies above it.
    Raises ValueError when a request cannot be met.
    """
    order = {node: i for i, node in enumerate(tree.traverse_preorder())}
    size: dict[treeswift.Node, int] = {}
    for node in tree.traverse_postorder():
        size[node] = (
            1 if node.is_leaf() else sum(size[child] for child in node.children)
        )

    has_heldout = {node: False for node in order}
    stem_parents: set[treeswift.Node] = set()
    groups: list[HeldOutGroup] = []

    def take(node, group_type, rank, taxon):
        members = tuple(leaf.label for leaf in node.traverse_leaves())
        groups.append(
            HeldOutGroup(
                f"{group_type}{sum(g.group_type == group_type for g in groups) + 1}",
                group_type,
                rank,
                taxon,
                members,
            )
        )
        stem_parents.add(node.parent)
        ancestor = node.parent
        while ancestor is not None and not has_heldout[ancestor]:
            has_heldout[ancestor] = True
            ancestor = ancestor.parent

    for rank in sorted(rank_counts):
        taxa = monophyletic_taxa(tree, lineages, rank)
        pool = sorted(
            (
                (node, taxon)
                for taxon, node in taxa.items()
                if min_leaves <= size[node] <= max_leaves
            ),
            key=lambda item: order[item[0]],
        )
        for _ in range(rank_counts[rank]):
            for i in rng.permutation(len(pool)):
                node, taxon = pool[i]
                if _is_free(node, has_heldout, stem_parents):
                    take(node, "clade", rank, taxon)
                    break
            else:
                raise ValueError(f"no eligible clade left at rank {rank}")

    leaves = [node for node in order if node.is_leaf()]
    for _ in range(n_singletons):
        for i in rng.permutation(len(leaves)):
            if _is_free(leaves[i], has_heldout, stem_parents):
                take(leaves[i], "singleton", None, leaves[i].label)
                break
        else:
            raise ValueError("no eligible singleton left")
    return groups


def sample_queries(
    groups: list[HeldOutGroup], rng: np.random.Generator, per_clade: int = 25
) -> list[str]:
    """Up to `per_clade` random members of each clade plus every singleton."""
    queries: list[str] = []
    for group in groups:
        if group.group_type == "singleton" or len(group.members) <= per_clade:
            queries.extend(group.members)
        else:
            keep = sorted(rng.choice(len(group.members), per_clade, replace=False))
            queries.extend(group.members[i] for i in keep)
    return queries


def prune_tree(tree: treeswift.Tree, heldout: set[str]) -> treeswift.Tree:
    """Backbone without `heldout`: unary nodes suppressed with additive lengths.

    The root edge is cleared and a degree-2 root contracted, as the backbone
    comparison in phylobench.accuracy expects.
    """
    unknown = heldout - {leaf.label for leaf in tree.traverse_leaves()}
    if unknown:
        raise ValueError(f"held-out labels not in tree: {sorted(unknown)[:5]}")
    with warnings.catch_warnings():
        # Suppression may delete internal labels; the backbone is label-free.
        warnings.filterwarnings("ignore", message="Deleting a node with label")
        backbone = tree.extract_tree_without(heldout)
    backbone.root.edge_length = None
    contract_bifurcating_root(backbone.root)
    return backbone


def prune_alignment(seqs: dict[str, str], heldout: set[str]) -> dict[str, str]:
    """Rows of `seqs` without `heldout`; every column is kept."""
    unknown = heldout - seqs.keys()
    if unknown:
        raise ValueError(f"held-out labels not in alignment: {sorted(unknown)[:5]}")
    return {name: seq for name, seq in seqs.items() if name not in heldout}
