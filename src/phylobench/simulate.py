"""Seeded GTR+Gamma4 simulator producing the tiny test fixture.

Files must regenerate byte-identically on any platform, so every operation on
simulated floats is an elementwise IEEE + - * / or a comparison. No matmul,
reductions, exp, log, pow or sqrt (BLAS/libm results can differ).
"""

import numpy as np
import treeswift

# GG 13_8 SEPP reference numbers (q2-fragment-insertion tests/data/ref-raxml-info.txt).
RATES = (
    0.947786,
    2.643649,
    1.480482,
    0.786790,
    4.471837,
    1.000000,
)  # ac ag at cg ct gt
FREQS = (0.250954, 0.234726, 0.318990, 0.195330)
# Discrete Gamma(4) mean-rate categories for alpha=1 (closed form).
GAMMA4_ALPHA1 = (0.136953, 0.476751, 1.0, 2.386294)
SEED = 2026
N_SITES = 1500
N_LEAVES = 110
PRIMER_BLOCKS = (
    (7, "AGAGTTTGATCCTGGCTCAG"),
    (337, "ACTCCTACGGGAGGCAGCAG"),
    (514, "GTGCCAGCAGCCGCGGTAA"),
    (783, "AGGATTAGATACCCTGGTAGTCC"),
    (1062, "AGGTGCTGCATGGCTGTCG"),
)
_MIN_SITES = 1081
_STATES = "ACGT"
_SQUARINGS = 10
_TAYLOR_ORDER = 8


def random_tree(n_leaves, rng, min_branch=0.005, max_branch=0.055):
    """Tree with a trifurcating root; branch lengths on a 1e-6 grid."""
    min_units = round(min_branch * 1_000_000)
    span_units = round(max_branch * 1_000_000) - min_units

    def new_child(parent):
        length = (min_units + int(rng.integers(0, span_units))) / 1_000_000
        node = treeswift.Node(edge_length=length)
        parent.add_child(node)
        return node

    root = treeswift.Node()
    leaves = [new_child(root) for _ in range(3)]
    while len(leaves) < n_leaves:
        i = int(rng.integers(len(leaves)))
        leaves[i], leaves[-1] = leaves[-1], leaves[i]
        parent = leaves.pop()
        leaves.extend([new_child(parent), new_child(parent)])
    tree = treeswift.Tree()
    tree.root = root
    for i, leaf in enumerate(tree.traverse_preorder(internal=False), start=1):
        leaf.label = f"t{i}"
    return tree


def gtr_matrix(rates, freqs):
    """GTR rate matrix normalized to mean rate 1, as nested lists."""
    q = [[0.0] * 4 for _ in range(4)]
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for r, (i, j) in zip(rates, pairs, strict=True):
        q[i][j] = r * freqs[j]
        q[j][i] = r * freqs[i]
    for i in range(4):
        total = 0.0
        for j in range(4):
            if j != i:
                total += q[i][j]
        q[i][i] = -total
    mean_rate = 0.0
    for i in range(4):
        mean_rate += freqs[i] * -q[i][i]
    return [[x / mean_rate for x in row] for row in q]


def _mm(a, b):
    out = a[..., :, 0, None] * b[..., None, 0, :]
    for k in range(1, 4):
        out = out + a[..., :, k, None] * b[..., None, k, :]
    return out


def transition_matrices(lengths, q):
    """P(rate_k * t) per edge and Gamma category, shape (E, 4, 4, 4).

    Scaling-and-squaring Taylor series built only from + * /.
    """
    lengths = np.asarray(lengths, dtype=np.float64)
    scale = np.array(GAMMA4_ALPHA1)[None, :] * lengths[:, None] / 2**_SQUARINGS
    m = scale[:, :, None, None] * np.array(q)[None, None]
    ident = np.eye(4)
    result = ident + m / _TAYLOR_ORDER
    for n in range(_TAYLOR_ORDER - 1, 0, -1):
        result = ident + _mm(m, result) / n
    for _ in range(_SQUARINGS):
        result = _mm(result, result)
    return result


def _sample_states(probs, u):
    c0 = probs[..., 0]
    c01 = c0 + probs[..., 1]
    c012 = c01 + probs[..., 2]
    return (
        (u > c0).astype(np.int8)
        + (u > c01).astype(np.int8)
        + (u > c012).astype(np.int8)
    )


def simulate_alignment(tree, n_sites, rng):
    """Gap-free GTR+Gamma4 alignment; returns {leaf label: sequence}."""
    nodes = list(tree.traverse_preorder())
    position = {id(node): i for i, node in enumerate(nodes)}
    lengths = [0.0 if node.is_root() else node.edge_length for node in nodes]
    p = transition_matrices(lengths, gtr_matrix(RATES, FREQS))
    category = rng.integers(0, 4, n_sites)

    root_probs = np.broadcast_to(np.array(FREQS), (n_sites, 4))
    states = [_sample_states(root_probs, rng.random(n_sites))]
    for i, node in enumerate(nodes[1:], start=1):
        parent_states = states[position[id(node.parent)]]
        states.append(
            _sample_states(p[i][category, parent_states], rng.random(n_sites))
        )

    alphabet = np.array(list(_STATES))
    return {
        node.label: "".join(alphabet[states[i]])
        for i, node in enumerate(nodes)
        if node.is_leaf()
    }


def embed_primer_sites(seqs):
    """Overwrite the four primer blocks in every sequence."""
    out = {}
    for name, seq in seqs.items():
        if len(seq) < _MIN_SITES:
            raise ValueError(f"{name}: alignment shorter than {_MIN_SITES} sites")
        chars = list(seq)
        for start, primer in PRIMER_BLOCKS:
            chars[start : start + len(primer)] = primer
        out[name] = "".join(chars)
    return out


def choose_heldout(tree, clade_sizes, n_singletons, rng):
    """Pick held-out groups as [(group_id, group_type, member labels)].

    A group is eligible only if its stem-parent subtree is otherwise untouched,
    so every group maps to its own distinct backbone bipartition.
    """
    nodes = list(tree.traverse_preorder())
    leaf_count = {}
    for node in reversed(nodes):
        leaf_count[id(node)] = (
            1 if node.is_leaf() else sum(leaf_count[id(c)] for c in node.children)
        )
    requests = [("clade", k, f"clade{i}") for i, k in enumerate(clade_sizes, start=1)]
    requests += [("singleton", 1, f"singleton{i}") for i in range(1, n_singletons + 1)]

    groups = []
    heldout = set()
    stem_parents = set()
    for group_type, size, group_id in requests:
        has_heldout = {}
        for node in reversed(nodes):
            has_heldout[id(node)] = (node.is_leaf() and node.label in heldout) or any(
                has_heldout[id(c)] for c in node.children
            )
        under_stem_parent = {id(nodes[0]): False}
        for node in nodes[1:]:
            parent = node.parent
            under_stem_parent[id(node)] = (
                under_stem_parent[id(parent)] or id(parent) in stem_parents
            )

        candidates = [
            node
            for node in nodes[1:]
            if leaf_count[id(node)] == size
            and (node.is_leaf() == (group_type == "singleton"))
            and not node.parent.is_root()
            and not has_heldout[id(node.parent)]
            and not under_stem_parent[id(node)]
        ]
        if not candidates:
            raise ValueError(f"no eligible {group_type} of size {size}")
        chosen = candidates[int(rng.integers(len(candidates)))]
        members = [leaf.label for leaf in chosen.traverse_preorder(internal=False)]
        groups.append((group_id, group_type, members))
        heldout.update(members)
        stem_parents.add(id(chosen.parent))
    return groups


def format_newick(tree):
    """One-line Newick: leaf labels only, 6-decimal lengths, no root length."""
    text = {}
    for node in tree.traverse_postorder():
        if node.is_leaf():
            body = node.label
        else:
            body = "(" + ",".join(text.pop(id(c)) for c in node.children) + ")"
        text[id(node)] = body if node.is_root() else f"{body}:{node.edge_length:.6f}"
    return text[id(tree.root)] + ";\n"


def format_fasta(seqs):
    return "".join(f">{name}\n{seq}\n" for name, seq in seqs.items())


def format_raxml_info(n_patterns):
    # Layout mirrors the q2-fragment-insertion/SEPP test file because pplacer only
    # parses this layout for RAxML 8; GG2 2022.10 layout is UNVERIFIED (T18 confirms).
    freqs = " ".join(f"{f:.6f}" for f in FREQS)
    rates = " ".join(f"{r:.6f}" for r in RATES)
    call = (
        "raxmlHPC-PTHREADS -s aligned-dna-sequences.fasta -m GTRCAT -n info "
        "-F -g tree.nwk -T 16 -p 3"
    )
    lines = [
        "This is RAxML version 8.1.24 released by Alexandros Stamatakis "
        "on June 29 2015.",
        "",
        f"Alignment has {n_patterns} distinct alignment patterns",
        "",
        "Partition: 0",
        f"Alignment Patterns: {n_patterns}",
        "Name: No Name Provided",
        "DataType: DNA",
        "Substitution Matrix: GTR",
        "",
        "RAxML was called as follows:",
        "",
        call,
        "",
        f"Base frequencies: {freqs} ",
        f"alpha[0]: 1.000000 rates[0] ac ag at cg ct gt: {rates}",
        "",
    ]
    return "\n".join(lines)


def format_truth(groups):
    rows = ["query\tgroup\tgroup_type"]
    for group_id, group_type, members in groups:
        rows += [f"{m}\t{group_id}\t{group_type}" for m in members]
    return "\n".join(rows) + "\n"


def tiny_fixture(seed=SEED):
    """Fixture files as {filename: text}.

    queries.fasta equals queries.aligned.fasta (gap-free), so the T7 projection
    of unaligned queries onto the backbone is not exercised by this fixture.
    """
    rng = np.random.default_rng(seed)
    reference = random_tree(N_LEAVES, rng)
    seqs = embed_primer_sites(simulate_alignment(reference, N_SITES, rng))
    groups = choose_heldout(reference, (3, 3), 4, rng)
    query_names = [m for _, _, members in groups for m in members]
    held = set(query_names)
    backbone = reference.extract_tree_without(held)
    backbone_seqs = {
        leaf.label: seqs[leaf.label]
        for leaf in reference.traverse_preorder(internal=False)
        if leaf.label not in held
    }
    query_seqs = {name: seqs[name] for name in query_names}
    n_patterns = len(set(zip(*backbone_seqs.values(), strict=True)))
    return {
        "reference_tree.nwk": format_newick(reference),
        "tree.nwk": format_newick(backbone),
        "aligned-dna-sequences.fasta": format_fasta(backbone_seqs),
        "queries.fasta": format_fasta(query_seqs),
        "queries.aligned.fasta": format_fasta(query_seqs),
        "raxml-info.txt": format_raxml_info(n_patterns),
        "truth.tsv": format_truth(groups),
    }
