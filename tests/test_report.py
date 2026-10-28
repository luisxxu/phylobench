import hashlib
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from phylobench.report import (
    bin_table,
    build_blocks,
    ci_reading,
    concordance_table,
    cost_table,
    distribution_table,
    extract_blocks,
    hardware_table,
    license_table,
    markdown_table,
    oracle_table,
    pass_count_table,
    paste_blocks,
    runs_overview,
    verdict_table,
)
from phylobench.summaries import (
    bin_summary,
    cell_distribution,
    end_to_end,
    oracle_gap,
    verdict,
)

NAN = math.nan
ROOT = Path(__file__).parents[1]
SAMPLE = ROOT / "tests" / "fixtures" / "sample_results"
SCRIPT = ROOT / "workflow" / "scripts" / "report_tables.py"
REAL_TABLES = ("runs", "accuracy", "concordance")
RESULTS = Path(os.environ.get("PHYLOBENCH_RESULTS", ROOT / "results"))


def read(name, **kwargs):
    return pd.read_csv(SAMPLE / name, sep="\t", **kwargs)


@pytest.fixture(scope="module")
def sample_accuracy():
    return read("accuracy.tsv")


@pytest.fixture(scope="module")
def sample_runs():
    return read("runs.tsv")


def acc(run_id, edge_dist, bins=None):
    n = len(edge_dist)
    return pd.DataFrame(
        {
            "run_id": run_id,
            "query": [f"q{i}" for i in range(n)],
            "sim": 0.9,
            "bin": bins or ["90-95"] * n,
            "edge_dist": edge_dist,
            "exact": [float(d == 0) for d in edge_dist],
            "delta_err": NAN,
            "label": "accuracy",
        }
    )


def cat(*frames):
    return pd.concat(frames, ignore_index=True)


def frame(**cols):
    return pd.DataFrame(cols)


def test_markdown_table_header_separator_and_escaping():
    text = markdown_table(frame(arm=["a|b"], n=[3]), index=False)
    lines = text.splitlines()
    assert lines[0] == "| arm | n |"
    assert lines[1] == "|---|---|"
    assert lines[2] == "| a\\|b | 3 |"


def test_markdown_table_keeps_index_as_first_column():
    df = frame(x=["1"], y=["2"], arm=["sepp-prod"]).set_index("arm")
    assert markdown_table(df).splitlines()[0] == "| arm | x | y |"


@pytest.mark.parametrize(
    ("lo", "hi", "reading"),
    [
        (0.1, 0.5, "worse"),
        (-0.5, -0.1, "better"),
        (-0.5, 0.5, "spans 0"),
        (0.0, 0.5, "spans 0"),
        (-0.5, 0.0, "spans 0"),
    ],
)
def test_ci_reading(lo, hi, reading):
    assert ci_reading(lo, hi) == reading


def two_arm_verdict():
    ref = acc("gg2/B10k/full-length/raw/sepp-prod", [0.0, 1.0, 2.0, 3.0])
    good = acc("gg2/B10k/full-length/hmmalign/pplacer", [0.0, 1.0, 2.0, 3.0])
    bad = acc("gg2/B10k/full-length/hmmalign/apples2", [9.0, 9.0, 9.0, 9.0])
    other_region = acc("gg2/B10k/V4-150/raw/sepp-prod", [0.0, 0.0])
    other_tool = acc("gg2/B10k/V4-150/hmmalign/pplacer", [0.0, 0.0])
    return verdict(cat(ref, good, bad, other_region, other_tool))


def test_verdict_table_cells_and_n_row():
    table = verdict_table(two_arm_verdict(), "gg2", "B10k")
    assert list(table.columns) == ["full-length", "V4-150"]
    assert table.loc["pplacer", "full-length"] == "pass +0.00 [+0.00, +0.00] spans 0"
    assert table.loc["apples2", "full-length"].startswith("fail +7.50")
    assert table.loc["apples2", "full-length"].endswith("worse")
    assert table.loc["n", "full-length"] == "4"
    assert table.loc["n", "V4-150"] == "2"


def test_verdict_table_blank_where_arm_absent_by_design():
    table = verdict_table(two_arm_verdict(), "gg2", "B10k")
    assert table.loc["apples2", "V4-150"] == ""


def test_verdict_table_states_without_statistics():
    ref = acc("gg2/Bfull/V4-150/raw/sepp-prod", [0.0, 1.0])
    unplaced = acc("gg2/Bfull/V4-150/hmmalign/pplacer", [NAN, NAN])
    table = verdict_table(verdict(cat(ref, unplaced)), "gg2", "Bfull")
    assert table.loc["pplacer", "V4-150"] == "unplaced"


def test_pass_count_table_counts_passes_and_large_regions():
    v = pd.DataFrame(
        {
            "dataset": "gg2",
            "backbone": "B10k",
            "region": ["a", "b", "c"],
            "arm": "pplacer",
            "n": [1000, 1000, 50],
            "verdict": ["pass", "fail", "pass"],
        }
    )
    table = pass_count_table(v)
    assert table.loc["pplacer", "gg2 B10k"] == "2/3 (1/2)"


def test_cell_distribution_median_fraction_at_zero_and_unplaced():
    df = acc("gg2/B10k/V4-150/hmmalign/pplacer", [0.0, 0.0, 3.0, NAN])
    row = cell_distribution(df).iloc[0]
    assert row["n"] == 4
    assert row["n_unplaced"] == 1
    assert row["median"] == 0.0
    assert row["mean"] == pytest.approx(1.0)
    assert row["frac_zero"] == pytest.approx(2 / 3)


def test_distribution_table_marks_all_unplaced():
    df = cat(
        acc("gg2/Bfull/V4-150/raw/sepp-prod", [0.0, 1.0]),
        acc("gg2/Bfull/V4-150/hmmalign/pplacer", [NAN, NAN]),
    )
    table = distribution_table(cell_distribution(df), "gg2", ["V4-150"])
    row = table[table["arm"] == "pplacer"].iloc[0]
    assert row["mean"] == "unplaced"
    assert row["unplaced"] == "2"


def test_bin_table_marks_small_bins_and_adds_n_row():
    big = acc(
        "gg2/B10k/V4-150/hmmalign/pplacer",
        [1.0] * 120 + [2.0] * 10,
        ["75-80"] * 120 + ["80-85"] * 10,
    )
    ref = acc(
        "gg2/B10k/V4-150/raw/sepp-prod",
        [0.0] * 120 + [0.0] * 10,
        ["75-80"] * 120 + ["80-85"] * 10,
    )
    bins = bin_summary(cat(big, ref))
    table = bin_table(bins, "gg2", "B10k", "V4-150")
    assert list(table.columns) == ["75-80", "80-85"]
    assert table.loc["pplacer", "75-80"] == "1.00"
    assert table.loc["pplacer", "80-85"] == "2.00†"
    assert table.loc["n", "80-85"] == "10"


def test_runs_overview_counts_status_and_memory_source():
    runs = frame(
        phase=["place", "place", "convert"],
        status=["ok", "DNF", "ok"],
        mem_source=["cgroup", "cgroup", "psutil"],
    )
    table = runs_overview(runs).set_index("phase")
    assert table.loc["place", "rows"] == "2"
    assert table.loc["place", "DNF"] == "1"
    assert table.loc["place", "cgroup"] == "2"
    assert table.loc["convert", "psutil"] == "1"
    assert table.loc["convert", "cgroup"] == "0"


def test_hardware_table_counts_cells_per_arm_backbone_and_cpu():
    e2e = frame(
        backbone=["B10k", "Bfull", "Bfull"],
        arm=["sepp-prod", "sepp-prod", "bscampp"],
        cpu_model=[
            "AMD EPYC 7302 16-Core Processor",
            "Intel(R) Xeon(R) Gold 6242 CPU @ 2.80GHz",
            "AMD EPYC 7302 16-Core Processor",
        ],
    )
    table = hardware_table(e2e)
    assert table.loc["sepp-prod", "B10k EPYC 7302"] == "1"
    assert table.loc["sepp-prod", "Bfull Xeon 6242"] == "1"
    assert table.loc["bscampp", "B10k EPYC 7302"] == "0"


def test_oracle_table_cells(sample_accuracy):
    table = oracle_table(oracle_gap(sample_accuracy), "gg2", "B10k")
    assert set(table.columns) <= {"full-length", "V4-150"}
    assert "[" in table.iloc[0, 0]


def test_cost_table_keeps_dnf_convert_failed_and_hardware(sample_runs):
    e2e = end_to_end(sample_runs)
    table = cost_table(e2e, "gg2", "Bfull", "full-length")
    assert (table["arm"] == "sepp-prod").any()
    flagged = table[table["e2e status"].str.contains("convert")]
    assert len(flagged) >= 1
    assert flagged["e2e vs sepp-prod"].eq("").all()
    assert {"arm", "peak mem MB", "mem source", "cpu"} <= set(table.columns)


def test_cost_table_ratio_against_reference():
    e2e = pd.DataFrame(
        {
            "dataset": "gg2",
            "backbone": "B10k",
            "region": "full-length",
            "aln": ["raw", "hmmalign"],
            "arm": ["sepp-prod", "pplacer"],
            "query_n": 1000.0,
            "place_s": [100.0, 25.0],
            "e2e_s": [100.0, 50.0],
            "peak_mem_mb": [1000.0, 250.0],
            "mem_source": "cgroup",
            "gpu_peak_mb": NAN,
            "cpu_model": "AMD EPYC 7302 16-Core Processor",
            "place_status": "ok",
            "convert_failed": False,
            "e2e_status": "ok",
        }
    )
    table = cost_table(e2e, "gg2", "B10k", "full-length").set_index("arm")
    assert table.loc["pplacer", "e2e vs sepp-prod"] == "0.50"
    assert table.loc["pplacer", "mem vs sepp-prod"] == "0.25"
    assert table.loc["pplacer", "cpu"] == "EPYC 7302"


def test_concordance_table_counts_nan_as_unplaced():
    df = pd.DataFrame(
        {
            "run_id": "emp-10k/B10k/V4-150/hmmalign/pplacer",
            "edge_dist": [0.0, 0.0, 4.0, NAN],
        }
    )
    row = concordance_table(df).iloc[0]
    assert row["dataset"] == "emp-10k"
    assert row["queries"] == "4"
    assert row["unplaced"] == "1"
    assert row["fraction at 0"] == "0.67"


LICENSES = """# Licenses

| Arm / tool | Version | Own license | Engine(s) executed | Copyleft flag | Source |
|---|---|---|---|---|---|
| `sepp-prod` | 2026.7.0 | BSD-3 | SEPP (GPL-3) | GPL-3 engines | https://x |
| HMMER | 3.4 | BSD-3 | itself | none flagged | https://y |

## Datasets
"""


def test_license_table_selects_columns_and_stops_at_blank_line():
    table = license_table(LICENSES)
    assert list(table.columns) == [
        "Arm / tool",
        "Version",
        "Own license",
        "Engine(s) executed",
        "Copyleft flag",
    ]
    assert len(table) == 2


def test_blocks_roundtrip_through_paste():
    blocks = {"a": "<!-- table:a -->\n**Table 1.** x\n\nNEW\n<!-- /table -->"}
    old = "intro\n<!-- table:a -->\nOLD\n<!-- /table -->\nouttro\n"
    pasted = paste_blocks(old, blocks)
    assert pasted == f"intro\n{blocks['a']}\nouttro\n"
    assert extract_blocks(pasted) == blocks


def test_paste_blocks_fails_on_missing_marker():
    with pytest.raises(ValueError, match="table:a"):
        paste_blocks("no markers", {"a": "<!-- table:a -->\nx\n<!-- /table -->"})


def test_build_blocks_on_sample_without_concordance(sample_runs, sample_accuracy):
    blocks = build_blocks(sample_runs, sample_accuracy, None, LICENSES)
    assert "verdict-gg2-B10k" in blocks
    assert not any(key.startswith("concordance") for key in blocks)
    assert all(re.search(r"\*\*Table \d+\.", block) for block in blocks.values())


def test_script_exits_zero_on_sample_tables_without_concordance():
    result = run_report_script(SAMPLE)
    assert result.returncode == 0, result.stderr[-3000:]
    assert "<!-- table:verdict-gg2-B10k -->" in result.stdout


def run_report_script(results):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--results",
            str(results),
            "--licenses",
            str(ROOT / "LICENSES.md"),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )


def header_hashes(text):
    return dict(
        (name, digest)
        for digest, name in re.findall(r"`([0-9a-f]{64})`\s+`(\w+\.tsv)`", text)
    )


@pytest.mark.skipif(
    not all((RESULTS / f"{n}.tsv").exists() for n in REAL_TABLES),
    reason="real result tables absent (set PHYLOBENCH_RESULTS); CI never has them",
)
def test_results_md_is_current_with_real_tables():
    text = (ROOT / "results.md").read_text()
    recorded = header_hashes(text)
    for name in REAL_TABLES:
        digest = hashlib.sha256((RESULTS / f"{name}.tsv").read_bytes()).hexdigest()
        assert recorded[f"{name}.tsv"] == digest, (
            f"{name}.tsv differs from results.md header"
        )
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--results",
            str(RESULTS),
            "--licenses",
            str(ROOT / "LICENSES.md"),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert result.returncode == 0, result.stderr[-3000:]
    assert extract_blocks(text) == extract_blocks(result.stdout)
