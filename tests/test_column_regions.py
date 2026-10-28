import subprocess
import sys
from pathlib import Path

import pytest
from test_regions import synthetic

from phylobench.regions import (
    COLUMN_REGIONS,
    COVERAGE_TOLERANCE,
    PRIMERS,
    REGIONS,
    ColumnWindow,
    column_window,
    cut_column_window,
    exclusion_counts,
    read_fasta,
    split_spans,
    ungap,
)

ROOT = Path(__file__).parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "tiny"
ALIGNED = dict(read_fasta(FIXTURE / "queries.aligned.fasta"))
QUERIES = dict(read_fasta(FIXTURE / "queries.fasta"))


def aligned_row(region, left=3, insert="A" * 200, width=400):
    """A synthetic aligned row: the primed sequence placed at column `left`."""
    seq, _ = synthetic(region, insert=insert)
    row = "-" * left + seq
    return row + "-" * (width - len(row))


def test_column_regions_name_every_primed_region_and_leave_the_originals():
    assert COLUMN_REGIONS == tuple(f"{r}c" for r in PRIMERS)
    assert not set(COLUMN_REGIONS) & set(REGIONS)
    assert len(COLUMN_REGIONS) == 8


def test_window_spans_the_amplicon_columns_between_primers():
    row = aligned_row("V3-V4", left=3)
    forward = len(PRIMERS["V3-V4"][0])
    window = column_window({"a": row}, "V3-V4c")
    assert window == ColumnWindow(3 + 3 + forward, 3 + 3 + forward + 200, 200)


def test_window_uses_the_median_over_rows_with_primer_hits():
    rows = {
        f"r{i}": aligned_row("V3-V4", left=left) for i, left in enumerate([2, 3, 9])
    }
    window = column_window(rows, "V3-V4c")
    assert window.start == column_window({"a": rows["r1"]}, "V3-V4c").start


def test_window_ignores_rows_without_primers():
    rows = {"hit": aligned_row("V3-V4"), "no_primers": "ACGT" * 100}
    assert column_window(rows, "V3-V4c") == column_window(
        {"hit": rows["hit"]}, "V3-V4c"
    )


def test_window_requires_at_least_one_primer_hit():
    with pytest.raises(ValueError, match="V3-V4c"):
        column_window({"x": "ACGT" * 100}, "V3-V4c")


def test_window_length_is_the_median_residue_count_inside_the_window():
    gapped = aligned_row("V3-V4", insert="A" * 100 + "-" * 5 + "A" * 95)
    rows = {
        "full": aligned_row("V3-V4"),
        "gapped": gapped,
        "full2": aligned_row("V3-V4"),
    }
    assert column_window(rows, "V3-V4c").length == 200


def test_v4_150c_window_is_150_columns_after_515f():
    row = aligned_row("V4-150", insert="ACGT" * 60)
    window = column_window({"a": row}, "V4-150c")
    assert (window.end - window.start, window.length) == (150, 150)
    assert row[window.start : window.end] == ("ACGT" * 60)[:150]


def window_of(length, start=10):
    return ColumnWindow(start, start + length, length)


def row_with_coverage(window, residues, width=200):
    """Residues fill the window from its left edge; the rest of the window is gaps."""
    cells = "A" * residues + "-" * (window.end - window.start - residues)
    return "-" * window.start + cells + "-" * (width - window.end)


def test_tolerance_is_a_fixed_constant():
    assert COVERAGE_TOLERANCE == 0.9


def test_query_at_exactly_the_tolerance_is_kept():
    window = window_of(100)
    kept, excluded = cut_column_window({"q": row_with_coverage(window, 90)}, window)
    assert list(kept) == ["q"] and excluded == {}
    assert len(kept["q"]) == 90


def test_query_just_below_the_tolerance_is_excluded():
    window = window_of(100)
    kept, excluded = cut_column_window({"q": row_with_coverage(window, 89)}, window)
    assert kept == {}
    assert exclusion_counts(excluded) == {"window coverage below tolerance": 1}


def test_query_without_primers_that_covers_the_window_is_kept():
    window = window_of(100)
    row = row_with_coverage(window, 100)
    assert split_spans([("q", ungap(row))], "V3-V4")[0] == {}
    kept, excluded = cut_column_window({"q": row}, window)
    assert list(kept) == ["q"] and excluded == {}


def test_cut_ignores_residues_outside_the_window_and_uppercases():
    window = window_of(10, start=5)
    row = "AAAAA" + "acgtacgtac" + "GGGGG"
    kept, _ = cut_column_window({"q": row}, window)
    assert kept["q"] == "ACGTACGTAC"


def test_exclusion_counts_cover_every_dropped_query():
    window = window_of(100)
    rows = {
        "ok": row_with_coverage(window, 100),
        "low1": row_with_coverage(window, 10),
        "low2": row_with_coverage(window, 0),
    }
    kept, excluded = cut_column_window(rows, window)
    assert list(kept) == ["ok"]
    assert exclusion_counts(excluded) == {"window coverage below tolerance": 2}


def test_original_regions_are_unchanged():
    assert REGIONS == (
        "full-length", "V4-150", "V3-V4", "V1-V2", "V1-V3", "V4-V5", "V5-V6",
        "V6-V8", "V7-V9",
    )  # fmt: skip
    spans, exclusions = split_spans(QUERIES.items(), "V4-150")
    assert exclusions == {} and len(spans) == len(QUERIES)


@pytest.mark.parametrize("region", COLUMN_REGIONS)
def test_fixture_queries_all_survive_every_column_region(region):
    window = column_window(ALIGNED, region)
    kept, excluded = cut_column_window(ALIGNED, window)
    assert excluded == {} and set(kept) == set(ALIGNED)
    for name, fragment in kept.items():
        assert fragment in QUERIES[name].upper()


def run_script(tmp_path, region, aligned=True):
    out = {n: tmp_path / n for n in ("q.fasta", "excl.tsv", "counts.tsv")}
    command = [
        sys.executable, str(ROOT / "workflow/scripts/extract_regions.py"),
        "--region", region, "--queries", str(FIXTURE / "queries.fasta"),
        "--out-fasta", str(out["q.fasta"]), "--out-exclusions", str(out["excl.tsv"]),
        "--out-counts", str(out["counts.tsv"]),
    ]  # fmt: skip
    if aligned:
        command += ["--aligned", str(FIXTURE / "queries.aligned.fasta")]
    subprocess.run(command, check=True, capture_output=True)
    return out


def test_script_writes_column_region_queries_and_counts(tmp_path):
    out = run_script(tmp_path, "V4-150c")
    written = dict(read_fasta(out["q.fasta"]))
    assert list(written) == list(QUERIES)
    assert (
        out["counts.tsv"].read_text() == f"reason\tcount\nextracted\t{len(QUERIES)}\n"
    )


def test_script_original_region_ignores_aligned_and_matches_amplicon_span(tmp_path):
    out = run_script(tmp_path, "V4-150", aligned=False)
    spans, _ = split_spans(QUERIES.items(), "V4-150")
    written = dict(read_fasta(out["q.fasta"]))
    assert written == {n: QUERIES[n].upper()[s:e] for n, (s, e) in spans.items()}
