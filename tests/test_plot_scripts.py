import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
SAMPLE = ROOT / "tests" / "fixtures" / "sample_results"
SCRIPTS = ROOT / "workflow" / "scripts"
ACCURACY_FIGURES = [
    "accuracy_edge_dist.png",
    "accuracy_exact.png",
    "accuracy_delta_err.png",
    "verdict.png",
    "oracle_gap.png",
]
COST_FIGURES = [
    "cost_vs_similarity.png",
    "cost_vs_backbone.png",
    "cost_vs_query_count.png",
    "end_to_end.png",
    "gpu_memory.png",
]


def run_script(script, outdir):
    args = ["--accuracy", str(SAMPLE / "accuracy.tsv"), "--outdir", str(outdir)]
    if script == "plot_cost.py":
        args += ["--runs", str(SAMPLE / "runs.tsv")]
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / script), *args],
        capture_output=True,
        text=True,
        env={**os.environ, "MPLBACKEND": "Agg", "PYTHONPATH": str(ROOT / "src")},
    )  # fmt: skip
    assert result.returncode == 0, result.stderr[-3000:]
    return result


@pytest.fixture(scope="module")
def figures(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("figures")
    run_script("plot_accuracy.py", outdir)
    run_script("plot_cost.py", outdir)
    return outdir


@pytest.mark.parametrize("name", ACCURACY_FIGURES + COST_FIGURES)
def test_every_figure_is_written(figures, name):
    assert (figures / name).stat().st_size > 10_000


def test_scripts_write_only_their_own_figures(tmp_path):
    run_script("plot_accuracy.py", tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(ACCURACY_FIGURES)


def test_scripts_create_the_outdir(tmp_path):
    run_script("plot_cost.py", tmp_path / "nested" / "figures")
    assert (tmp_path / "nested" / "figures" / "gpu_memory.png").exists()


@pytest.fixture(scope="module")
def built():
    sys.path.insert(0, str(SCRIPTS))
    import matplotlib

    matplotlib.use("Agg")
    import pandas as pd
    import plot_accuracy
    import plot_cost

    from phylobench.summaries import bin_summary, end_to_end

    runs = pd.read_csv(SAMPLE / "runs.tsv", sep="\t")
    accuracy = pd.read_csv(SAMPLE / "accuracy.tsv", sep="\t")
    cells = end_to_end(runs)
    summary = bin_summary(accuracy)
    both = ["B10k", "Bfull"]
    figures = {
        "edge": plot_accuracy.metric_figure(summary, "mean_edge_dist", "e", both),
        "exact": plot_accuracy.metric_figure(summary, "exact_frac", "x", both),
        "delta": plot_accuracy.metric_figure(summary, "mean_delta_err", "d", ["B10k"]),
        "verdict": plot_accuracy.verdict_figure(accuracy),
        "oracle": plot_accuracy.oracle_gap_figure(accuracy),
        "similarity": plot_cost.cost_vs_similarity(
            runs, cells, accuracy[["run_id", "sim"]]
        ),
        "query_count": plot_cost.cost_vs_query_count(runs, cells),
        "e2e": plot_cost.end_to_end_figure(runs, cells),
        "gpu": plot_cost.gpu_memory(runs, cells),
    }
    yield figures, cells
    import matplotlib.pyplot as plt

    plt.close("all")


def visible_titles(fig):
    return [ax.get_title() for ax in fig.axes if ax.get_visible()]


@pytest.mark.parametrize(
    "name",
    [
        "edge",
        "exact",
        "delta",
        "verdict",
        "oracle",
        "similarity",
        "query_count",
        "e2e",
        "gpu",
    ],
)
def test_no_panel_mixes_backbones(built, name):
    titles = visible_titles(built[0][name])
    assert titles
    for title in titles:
        assert sum(b in title.split() for b in ("B10k", "Bfull")) == 1, title


def test_delta_error_figure_is_b10k_only(built):
    assert all("B10k" in t for t in visible_titles(built[0]["delta"]))


def test_unplaced_cell_is_annotated_with_placed_count(built):
    texts = [t.get_text() for ax in built[0]["edge"].axes for t in ax.texts]
    assert "pplacer: 0/24 placed" in texts
    verdict_texts = [t.get_text() for ax in built[0]["verdict"].axes for t in ax.texts]
    assert "0/24 placed" in verdict_texts


def dnf_markers(fig):
    return sum(
        len(c.get_offsets())
        for ax in fig.axes
        for c in ax.collections
        if c.get_gid() == "dnf"
    )


def test_end_to_end_marks_every_dnf_cell(built):
    figures, cells = built
    plotted = cells.dropna(subset=["place_s", "e2e_s"])
    assert dnf_markers(figures["e2e"]) == (plotted["e2e_status"] == "DNF").sum() > 0


def test_cost_figures_mark_placement_dnf(built):
    figures, cells = built
    failed = (cells["place_status"] == "DNF") & cells["dataset"].eq("gg2")
    assert failed.sum() == 1
    assert dnf_markers(figures["similarity"]) == 2 * failed.sum()


def test_cost_figures_carry_hardware_caveat(built):
    for name in ("similarity", "query_count", "e2e", "gpu"):
        texts = [t.get_text() for t in built[0][name].texts]
        assert any("Hardware:" in t and "sepp-prod" in t for t in texts), name


def test_similarity_figure_captions_the_per_cell_proxy(built):
    texts = [t.get_text() for t in built[0]["similarity"].texts]
    assert any("per cell, not per similarity bin" in t for t in texts)
