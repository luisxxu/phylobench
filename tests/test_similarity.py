from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phylobench.similarity import bin_report, similarity, similarity_bin

FIXTURE = Path(__file__).parent / "fixtures" / "tiny"


def read_fasta(path):
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line
    return seqs


def brute_force(query, refs, window):
    start, end = window
    columns = [i for i in range(start, end) if query[i] not in "-."]
    best = max(
        sum(
            1
            for i in columns
            if query[i].upper().replace("U", "T") in "ACGT"
            and query[i].upper().replace("U", "T") == ref[i].upper().replace("U", "T")
        )
        for ref in refs
    )
    return best, len(columns)


def row(df, query):
    r = df.set_index("query").loc[query]
    return int(r.matches), int(r.n_columns), float(r.sim)


def test_hand_built_rules():
    refs = {"r1": "ACGTACGT", "r2": "ACG-ACGT", "r3": "NNNNNNNN"}
    queries = {"q": "AC-TAcgN", "gaponly": "ACG.----"}
    df = similarity(queries, refs)
    assert list(df.columns) == ["query", "matches", "n_columns", "sim"]
    # q: non-gap columns 0,1,3,4,5,6,7 (N counts, never matches); r1 matches 6
    assert row(df, "q") == (6, 7, 6 / 7)
    assert row(df, "gaponly") == (3, 3, 1.0)


def test_ref_gap_and_ambiguous_are_mismatches_and_max_over_refs():
    refs = {"r1": "A-NT", "r2": "ACGA", "r3": "TCGT"}
    df = similarity({"q": "ACGT"}, refs)
    assert row(df, "q") == (3, 4, 0.75)
    df = similarity({"q": "ACGT"}, {"r1": "A-NT"})
    assert row(df, "q") == (2, 4, 0.5)


def test_u_is_t_and_case_insensitive():
    df = similarity({"q": "acgu"}, {"r": "ACGT"})
    assert row(df, "q") == (4, 4, 1.0)


def test_windows_differ_per_query_in_one_chunk():
    refs = {"r": "AAAACCCC"}
    queries = {"a": "AAAAGGGG", "b": "TTTTCCCC"}
    df = similarity(queries, refs, windows={"a": (0, 4), "b": (4, 8)}, chunk=2)
    assert row(df, "a") == (4, 4, 1.0)
    assert row(df, "b") == (4, 4, 1.0)
    df = similarity(queries, refs, windows={"a": (2, 6), "b": (0, 6)}, chunk=2)
    assert row(df, "a") == (2, 4, 0.5)
    assert row(df, "b") == (2, 6, 2 / 6)


def test_empty_window_raises_naming_query():
    with pytest.raises(ValueError, match="qgap"):
        similarity({"ok": "ACGT", "qgap": "A--T"}, {"r": "ACGT"}, {"qgap": (1, 3)})
    with pytest.raises(ValueError, match="allgap"):
        similarity({"allgap": "----"}, {"r": "ACGT"})


@pytest.mark.parametrize("chunk", [1, 7, 1000])
def test_random_matches_brute_force_exactly(chunk):
    rng = np.random.default_rng(3)
    alphabet = np.array(list("ACGTACGTACGTACGT-.NRacgu"))
    n_cols = 60

    def seqs(prefix, n):
        return {f"{prefix}{i}": "".join(rng.choice(alphabet, n_cols)) for i in range(n)}

    refs, queries = seqs("r", 40), seqs("q", 30)
    windows = {}
    for name in list(queries)[:15]:
        start = int(rng.integers(0, 30))
        windows[name] = (start, start + int(rng.integers(20, 30)))
    df = similarity(queries, refs, windows, chunk=chunk)
    assert df["query"].tolist() == list(queries)
    for name, q in queries.items():
        expected = brute_force(q, list(refs.values()), windows.get(name, (0, n_cols)))
        assert row(df, name)[:2] == expected
        assert row(df, name)[2] == pytest.approx(expected[0] / expected[1])


@pytest.mark.parametrize(
    ("matches", "n_columns", "label"),
    [
        (999, 1000, ">=99.9"),
        (1000, 1000, ">=99.9"),
        (998, 1000, "99-99.9"),
        (990, 1000, "99-99.9"),
        (989, 1000, "97-99"),
        (970, 1000, "97-99"),
        (969, 1000, "95-97"),
        (950, 1000, "95-97"),
        (900, 1000, "90-95"),
        (850, 1000, "85-90"),
        (800, 1000, "80-85"),
        (750, 1000, "75-80"),
        (749, 1000, "<75"),
        (0, 5, "<75"),
    ],
)
def test_similarity_bin_boundaries(matches, n_columns, label):
    out = similarity_bin(np.array([matches]), np.array([n_columns]))
    assert out.tolist() == [label]


def test_similarity_bin_is_exact_at_boundary_with_odd_lengths():
    # 0.999 * 1501 = 1499.5 -> 1499 is below, 1500 is at or above
    out = similarity_bin(np.array([1499, 1500]), np.array([1501, 1501]))
    assert out.tolist() == ["99-99.9", ">=99.9"]


def test_bin_report_counts_flags_and_lengths():
    table = pd.DataFrame(
        {
            "query": [f"q{i}" for i in range(104)],
            "bin": [">=99.9"] * 100 + ["75-80"] * 4,
            "lineage": ["g1"] * 60 + ["g2"] * 40 + ["g3"] * 3 + ["g4"],
            "length": list(range(100, 200)) + [10, 20, 30, 40],
        }
    )
    composition, summary = bin_report(table)
    comp = {(r.bin, r.lineage): r.n for r in composition.itertuples()}
    assert comp == {
        (">=99.9", "g1"): 60,
        (">=99.9", "g2"): 40,
        ("75-80", "g3"): 3,
        ("75-80", "g4"): 1,
    }
    s = summary.set_index("bin")
    assert s.loc[">=99.9", "n"] == 100
    assert not s.loc[">=99.9", "flagged"]
    assert s.loc["75-80", "flagged"]
    assert (s.loc[">=99.9", "length_min"], s.loc[">=99.9", "length_max"]) == (100, 199)
    assert s.loc[">=99.9", "length_median"] == 149.5
    assert s.loc["75-80", "length_median"] == 25


def test_fixture_queries_against_backbone():
    refs = read_fasta(FIXTURE / "aligned-dna-sequences.fasta")
    queries = read_fasta(FIXTURE / "queries.aligned.fasta")
    df = similarity(queries, refs)
    assert len(df) == len(queries)
    assert ((df.sim > 0) & (df.sim <= 1)).all()
    labels = similarity_bin(df.matches.to_numpy(), df.n_columns.to_numpy())
    assert len(labels) == len(df)
    assert not (labels == "").any()
