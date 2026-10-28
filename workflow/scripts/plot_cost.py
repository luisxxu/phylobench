import argparse
from pathlib import Path

import pandas as pd
from plot_common import (
    arm_color,
    arm_labels,
    footnote,
    grid,
    hardware_note,
    legend,
    present_arms,
    save,
)

from phylobench.summaries import end_to_end, split_run_id

BACKBONES = ["B10k", "Bfull"]
LOG_LABELS = {
    "place_s": "placement wall time (s)",
    "peak_mem_mb": "peak host memory (MB)",
}
DNF_NOTE = (
    "x: DNF, plotted at the wall time or memory reached. "
    "Hollow markers: mem_source=psutil."
)
QUERY_COUNT_DATASETS = ["gg2", "emp-10k", "emp-100k"]
V4 = "V4-150"


def draw(ax, df, x, y, *, status="place_status", hollow_psutil=False, line=True):
    for arm, g in df.groupby("arm"):
        color = arm_color(arm)
        ok = g[g[status] == "ok"].dropna(subset=[x, y])
        if line and len(ok):
            median = ok.groupby(x)[y].median()
            ax.plot(median.index, median.values, color=color, lw=1)
        hollow = (
            ok["mem_source"].eq("psutil")
            if hollow_psutil
            else pd.Series(False, index=ok.index)
        )
        ax.scatter(ok.loc[~hollow, x], ok.loc[~hollow, y], color=color, s=12, zorder=3)
        ax.scatter(
            ok.loc[hollow, x],
            ok.loc[hollow, y],
            facecolors="white",
            edgecolors=color,
            s=12,
            zorder=3,
        )
        failed = g[g[status] == "DNF"].dropna(subset=[x, y])
        if len(failed):
            ax.scatter(
                failed[x], failed[y], marker="x", color=color, s=30, zorder=4, gid="dnf"
            )


def finish(fig, runs, df, *lines):
    legend(fig, present_arms(df["arm"]), arm_labels(runs))
    footnote(fig, *lines, DNF_NOTE, hardware_note(runs))


def cost_rows(fig_axes, df, x):
    for (metric, backbone), ax in fig_axes.items():
        draw(
            ax,
            df[df["backbone"] == backbone],
            x,
            metric,
            hollow_psutil=metric == "peak_mem_mb",
        )
        ax.set_yscale("log")
        ax.set_title(backbone, fontsize=9)
        if backbone == BACKBONES[0]:
            ax.set_ylabel(LOG_LABELS[metric], fontsize=8)


def metric_grid(df):
    backbones = [b for b in BACKBONES if b in set(df["backbone"])]
    return grid(["place_s", "peak_mem_mb"], backbones)


def cost_vs_similarity(runs, cells, accuracy):
    median_sim = (
        split_run_id(accuracy)
        .groupby(["dataset", "backbone", "region"])["sim"]
        .median()
        .rename("median_sim")
        .reset_index()
    )
    df = cells[cells["dataset"] == "gg2"].merge(
        median_sim, on=["dataset", "backbone", "region"]
    )
    fig, axes = metric_grid(df)
    cost_rows(axes, df, "median_sim")
    for ax in axes.values():
        ax.set_xlabel("median query similarity of the cell (fraction)", fontsize=8)
    finish(
        fig,
        runs,
        df,
        "Runtime and memory are measured per cell, not per similarity bin: "
        "each point is one gg2 cell (one region) at the region's median similarity.",
    )
    return fig


def cost_vs_backbone(runs, cells):
    datasets = [d for d in ["gg2", "emp-10k", "emp-100k"] if d in set(cells["dataset"])]
    df = cells[
        cells["dataset"].isin(datasets) & cells["backbone"].isin(BACKBONES)
    ].copy()
    arms = present_arms(df["arm"])
    offset = df["arm"].map({a: (i - len(arms) / 2) * 0.03 for i, a in enumerate(arms)})
    df["xpos"] = df["backbone"].map({b: i for i, b in enumerate(BACKBONES)}) + offset
    fig, axes = grid(["place_s", "peak_mem_mb"], datasets)
    for (metric, dataset), ax in axes.items():
        draw(
            ax,
            df[df["dataset"] == dataset],
            "xpos",
            metric,
            hollow_psutil=metric == "peak_mem_mb",
        )
        ax.set_yscale("log")
        ax.set_xticks(range(len(BACKBONES)), BACKBONES)
        ax.set_title(
            dataset if dataset != "gg2" else "gg2 (one point per region)", fontsize=9
        )
        if dataset == datasets[0]:
            ax.set_ylabel(LOG_LABELS[metric], fontsize=8)
    finish(
        fig, runs, df, "x: backbone size (B10k then Bfull); lines join per-arm medians."
    )
    return fig


def cost_vs_query_count(runs, cells):
    df = cells[(cells["region"] == V4) & cells["dataset"].isin(QUERY_COUNT_DATASETS)]
    fig, axes = metric_grid(df)
    cost_rows(axes, df, "query_n")
    for ax in axes.values():
        ax.set_xscale("log")
        ax.set_xlabel("queries placed (gg2 V4-150, emp-10k, emp-100k)", fontsize=8)
    finish(fig, runs, df, "All cells at V4-150; log-log.")
    return fig


def end_to_end_figure(runs, cells):
    df = cells.dropna(subset=["place_s", "e2e_s"])
    fig, axes = grid(["e2e"], [b for b in BACKBONES if b in set(df["backbone"])])
    for (_, backbone), ax in axes.items():
        cell = df[df["backbone"] == backbone]
        draw(ax, cell, "place_s", "e2e_s", status="e2e_status", line=False)
        lo, hi = cell["place_s"].min(), cell["e2e_s"].max()
        ax.plot([lo, hi], [lo, hi], color="0.5", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(backbone, fontsize=9)
        ax.set_xlabel("placement-only wall time (s)", fontsize=8)
        if backbone == BACKBONES[0]:
            ax.set_ylabel("end-to-end wall time (s)", fontsize=8)
    finish(
        fig,
        runs,
        df,
        "End to end = region hmmalign + place + convert (SEPP: its own phases); "
        "ref-prep and hmmbuild excluded. Dashed: y = x.",
        "x marks a failed convert (apples2 guppy tog, pplacer at gg2 Bfull V4-150) "
        "or a failed placement; accuracy is still scored where placed.",
    )
    return fig


def gpu_memory(runs, cells):
    df = cells[cells["arm"].isin(["dipper-m", "dipper-m-gpu"])].dropna(
        subset=["query_n"]
    )
    gpu = df[df["arm"] == "dipper-m-gpu"]
    backbones = [b for b in BACKBONES if b in set(df["backbone"])]
    fig, axes = grid(["gpu_peak_mb", "peak_mem_mb"], backbones)
    for (metric, backbone), ax in axes.items():
        source = gpu if metric == "gpu_peak_mb" else df
        draw(
            ax,
            source[source["backbone"] == backbone],
            "query_n",
            metric,
            hollow_psutil=metric == "peak_mem_mb",
        )
        ax.set_xscale("log")
        if metric != "gpu_peak_mb":
            ax.set_yscale("log")
        ax.set_title(backbone, fontsize=9)
        ax.set_xlabel("queries placed", fontsize=8)
        if backbone == backbones[0]:
            ax.set_ylabel(
                "GPU peak memory (MB, dipper-m-gpu)"
                if metric == "gpu_peak_mb"
                else "host peak memory (MB)",
                fontsize=8,
            )
    finish(
        fig,
        runs,
        df,
        "Top: GPU memory of dipper-m-gpu. "
        "Bottom: host memory of dipper-m and dipper-m-gpu.",
    )
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Runtime and memory figures from runs.tsv."
    )
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--accuracy", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    runs = pd.read_csv(args.runs, sep="\t")
    accuracy = pd.read_csv(args.accuracy, sep="\t", usecols=["run_id", "sim"])
    cells = end_to_end(runs)
    save(
        cost_vs_similarity(runs, cells, accuracy), args.outdir, "cost_vs_similarity.png"
    )
    save(cost_vs_backbone(runs, cells), args.outdir, "cost_vs_backbone.png")
    save(cost_vs_query_count(runs, cells), args.outdir, "cost_vs_query_count.png")
    save(end_to_end_figure(runs, cells), args.outdir, "end_to_end.png")
    save(gpu_memory(runs, cells), args.outdir, "gpu_memory.png")


if __name__ == "__main__":
    main()
