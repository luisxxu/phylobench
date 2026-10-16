import shutil
import subprocess
from pathlib import Path

import pytest

from phylobench.regions import (
    PRIMERS,
    REGIONS,
    TRIM,
    amplicon_span,
    check_ungaps_to,
    exclusion_counts,
    mask_inserts,
    primer_pattern,
    project,
    read_fasta,
    reverse_complement,
    split_spans,
    stockholm_lines,
    ungap,
)

FIXTURE = Path(__file__).parents[1] / "tests" / "fixtures" / "tiny"
QUERIES = dict(read_fasta(FIXTURE / "queries.fasta"))
BACKBONE = dict(read_fasta(FIXTURE / "aligned-dna-sequences.fasta"))
BACKBONE_WIDTH = 1500
PRIMED_REGIONS = [r for r in REGIONS if r != "full-length"]


def synthetic(region, insert="A" * 200):
    forward, reverse = PRIMERS[region]
    forward = forward.replace("Y", "C").replace("M", "A").replace("N", "G")
    forward = forward.replace("W", "A")
    reverse_site = reverse_complement(reverse).replace("N", "A").replace("B", "C")
    reverse_site = reverse_site.replace("W", "A").replace("V", "A").replace("H", "A")
    reverse_site = reverse_site.replace("D", "A")
    return f"GGG{forward}{insert}{reverse_site}TTT", insert


def test_reverse_complement_handles_iupac():
    assert reverse_complement("GTGYCAGCM") == "KGCTGRCAC"


def test_primer_pattern_matches_iupac_exactly():
    pattern = primer_pattern("GTGYCAGCM")
    assert pattern.search("AAGTGCCAGCAAA")
    assert pattern.search("AAGTGTCAGCCAA")
    assert not pattern.search("AAGTGGCAGCAAA")


def test_ungap_removes_gap_characters():
    assert ungap("A-C.G--T") == "ACGT"


def test_full_length_span_is_whole_sequence():
    assert amplicon_span("ACGTACGT", "full-length") == (0, 8)


@pytest.mark.parametrize("region", [r for r in PRIMED_REGIONS if r not in TRIM])
def test_span_excludes_primers(region):
    seq, insert = synthetic(region)
    start, end = amplicon_span(seq, region)
    assert seq[start:end] == insert


def test_v4_150_takes_first_150_after_forward_primer():
    seq, insert = synthetic("V4-150", insert="ACGT" * 60)
    start, end = amplicon_span(seq, "V4-150")
    assert seq[start:end] == insert[:150]
    assert end - start == TRIM["V4-150"]


def test_v4_150_requires_reverse_primer():
    seq, _ = synthetic("V4-150", insert="ACGT" * 60)
    forward_only = seq[: seq.index(reverse_complement(PRIMERS["V4-150"][1])[:6])]
    assert amplicon_span(forward_only, "V4-150") == "no reverse primer"


def test_v4_150_amplicon_shorter_than_trim_is_too_short():
    seq, _ = synthetic("V4-150", insert="ACGT" * 10)
    assert amplicon_span(seq, "V4-150") == "too short"


def test_missing_primers_return_reasons():
    assert amplicon_span("ACGT" * 100, "V3-V4") == "no forward primer"
    seq, _ = synthetic("V3-V4")
    no_reverse = seq[: seq.index(reverse_complement(PRIMERS["V3-V4"][1])[:6])]
    assert amplicon_span(no_reverse, "V3-V4") == "no reverse primer"


def test_reverse_primer_must_follow_forward_primer():
    seq, _ = synthetic("V3-V4")
    forward, reverse = PRIMERS["V3-V4"]
    reverse_site = seq[-3 - len(reverse) : -3]
    assert amplicon_span(reverse_site + "A" * 50, "V3-V4") == "no forward primer"


def test_mismatched_primer_is_not_matched():
    seq, _ = synthetic("V3-V4")
    broken = seq.replace("CCTACGGG", "CCTACGGC", 1)
    assert amplicon_span(broken, "V3-V4") == "no forward primer"


@pytest.mark.parametrize("region", REGIONS)
def test_every_region_extracts_for_every_fixture_query(region):
    spans, exclusions = split_spans(QUERIES.items(), region)
    assert exclusions == {}
    assert set(spans) == set(QUERIES)
    for name, (start, end) in spans.items():
        assert 0 <= start < end <= len(QUERIES[name])


@pytest.mark.parametrize("region", PRIMED_REGIONS)
def test_fixture_spans_exclude_primer_sites(region):
    spans, _ = split_spans(QUERIES.items(), region)
    forward, reverse = PRIMERS[region]
    for name, (start, end) in spans.items():
        seq = QUERIES[name]
        assert primer_pattern(forward).fullmatch(seq[start - len(forward) : start])
        if region != "V4-150":
            assert primer_pattern(reverse_complement(reverse)).match(seq, end)


def test_split_spans_reports_exclusion_reasons():
    seqs = {"ok": synthetic("V3-V4")[0], "bad": "ACGT" * 50}
    spans, exclusions = split_spans(seqs.items(), "V3-V4")
    assert list(spans) == ["ok"]
    assert exclusions == {"bad": "no forward primer"}


def test_exclusion_counts():
    exclusions = {"a": "too short", "b": "too short", "c": "no reverse primer"}
    assert exclusion_counts(exclusions) == {"too short": 2, "no reverse primer": 1}


def test_project_maps_ungapped_span_to_columns():
    aligned = "--AC-GT.TT-"
    assert project((0, 2), aligned) == (2, 4)
    assert project((1, 4), aligned) == (3, 7)
    assert project((2, 6), aligned) == (5, 10)


def test_project_identity_on_gapless_row():
    assert project((5, 9), "ACGTACGTACGT") == (5, 9)


def test_check_ungaps_to_accepts_matching_rows():
    check_ungaps_to({"a": "A-CG", "b": "-TT-"}, {"a": "ACG", "b": "TT"})


def test_check_ungaps_to_rejects_changed_sequence():
    with pytest.raises(ValueError, match="b"):
        check_ungaps_to({"a": "A-CG", "b": "-TA-"}, {"a": "ACG", "b": "TT"})


def test_check_ungaps_to_rejects_name_mismatch():
    with pytest.raises(ValueError, match="names"):
        check_ungaps_to({"a": "A-CG"}, {"a": "ACG", "b": "TT"})


def test_fixture_aligned_queries_ungap_to_full_length():
    aligned = dict(read_fasta(FIXTURE / "queries.aligned.fasta"))
    check_ungaps_to(aligned, QUERIES)


def test_stockholm_lines_have_all_x_rf_line():
    rows = [("a", "AC-GT"), ("b", "A--GT")]
    assert list(stockholm_lines(rows)) == [
        "# STOCKHOLM 1.0",
        "a AC-GT",
        "b A--GT",
        "#=GC RF xxxxx",
        "//",
    ]


def test_stockholm_lines_reject_ragged_rows():
    with pytest.raises(ValueError, match="b"):
        list(stockholm_lines([("a", "ACGT"), ("b", "AC")]))


def test_mask_inserts_deletes_lowercase_and_dots():
    assert mask_inserts("AC-acg..GT.t") == "AC-GT"


@pytest.mark.skipif(
    not (shutil.which("hmmbuild") and shutil.which("hmmalign")),
    reason="hmmbuild/hmmalign not on PATH",
)
@pytest.mark.parametrize("region", REGIONS)
def test_hmmalign_masked_rows_have_backbone_column_count(tmp_path, region):
    sto = tmp_path / "backbone.sto"
    sto.write_text("\n".join(stockholm_lines(BACKBONE.items())) + "\n")
    hmm = tmp_path / "backbone.hmm"
    subprocess.run(
        ["hmmbuild", "--dna", "--hand", "--cpu", "1", str(hmm), str(sto)],
        check=True,
        capture_output=True,
    )
    spans, _ = split_spans(QUERIES.items(), region)
    # one insertion and one deletion so the mask step is exercised
    names = list(spans)
    queries = {n: QUERIES[n][slice(*spans[n])] for n in names}
    mid = len(queries[names[0]]) // 2
    queries[names[0]] = queries[names[0]][:mid] + "GGGGG" + queries[names[0]][mid:]
    queries[names[1]] = queries[names[1]][: mid - 10] + queries[names[1]][mid:]
    fasta = tmp_path / "queries.fasta"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in queries.items()))
    afa = tmp_path / "out.afa"
    subprocess.run(
        [
            "hmmalign",
            "--dna",
            "--outformat",
            "afa",
            "-o",
            str(afa),
            str(hmm),
            str(fasta),
        ],
        check=True,
        capture_output=True,
    )
    aligned = {n: mask_inserts(s) for n, s in read_fasta(afa)}
    assert list(aligned) == names
    assert {len(s) for s in aligned.values()} == {BACKBONE_WIDTH}
