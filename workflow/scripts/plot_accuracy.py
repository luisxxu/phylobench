import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plot_common import (
    APPLES_CAVEAT,
    arm_color,
    footnote,
    hide_unused,
    legend,
    present_arms,
    present_regions,
    region_grid,
    save,
    side_legend,
)

from phylobench.summaries import (
    BIN_ORDER,
    MARGIN,
    SMALL_BIN_N,
    bin_summary,
    oracle_gap,
    verdict,
)

VERDICT_COLORS = {
    "pass": "tab:green",
    "fail": "tab:red",
    "unplaced": "tab:orange",
    "no-reference": "tab:gray",
}
VERDICT_ROWS = [("gg2", "B10k"), ("gg2", "Bfull"), ("sim10k", "B10k")]


def unplaced_notes(cell: pd.DataFrame) -> list[tuple[str, str]]:
    cell = cell.groupby("arm")[["n", "n_unplaced"]].sum()
    return [
        (arm, f"{arm}: {r.n - r.n_unplaced}/{r.n} placed")
        for arm, r in cell[cell["n_unplaced"] > 0].iterrows()
    ]


def metric_figure(summary, metric, ylabel, backbones):
    df = summary[
        (summary["dataset"] == "gg2")
        & (summary["aln"] != "oracle")
        & summary["backbone"].isin(backbones)
    ]
    backbones = [b for b in backbones if b in set(df["backbone"])]
    regions = present_regions(df["region"])
    fig, axes = region_grid(dict.fromkeys(backbones, regions))
    for (backbone, region), ax in axes.items():
        cell = df[(df["backbone"] == backbone) & (df["region"] == region)]
        for arm, g in cell.groupby("arm"):
            by_bin = g.set_index(g["bin"].astype(str)).reindex(BIN_ORDER)
            x = np.arange(len(BIN_ORDER))
            y = by_bin[metric].to_numpy(dtype=float)
            small = (by_bin["n"] < SMALL_BIN_N).to_numpy()
            color = arm_color(arm)
            ax.plot(x, y, color=color, lw=1)
            ax.plot(x[~small], y[~small], "o", color=color, ms=3)
            ax.plot(x[small], y[small], "o", color=color, mfc="white", ms=3)
        for i, (arm, note) in enumerate(unplaced_notes(cell)):
            ax.text(
                0.02,
                0.98 - 0.07 * i,
                note,
                transform=ax.transAxes,
                fontsize=6,
                va="top",
                color=arm_color(arm),
            )
        ax.set_title(f"{backbone} {region}", fontsize=9)
        ax.set_xticks(range(len(BIN_ORDER)), BIN_ORDER, rotation=90, fontsize=6)
        if ax.get_subplotspec().is_first_col():
            ax.set_ylabel(ylabel, fontsize=8)
    legend(fig, present_arms(df["arm"]), {})
    footnote(
        fig,
        f"x: query-to-backbone similarity bin (%). "
        f"Hollow markers: bin has n < {SMALL_BIN_N}.",
        APPLES_CAVEAT,
    )
    return fig


def interval_figure(df, keys, color_of, xlabel, notes, rule_lines=()):
    regions_by_key = {
        (dataset, backbone): present_regions(
            df.loc[(df["dataset"] == dataset) & (df["backbone"] == backbone), "region"]
        )
        for dataset, backbone in keys
    }
    regions_by_key = {k: r for k, r in regions_by_key.items() if r}
    arms = present_arms(df["arm"])[::-1]
    fig, axes = region_grid(regions_by_key, panel=(3.2, 0.3 * len(arms) + 1.2))
    used = set()
    for ((dataset, backbone), region), ax in axes.items():
        cell = df[
            (df["dataset"] == dataset)
            & (df["backbone"] == backbone)
            & (df["region"] == region)
        ]
        if cell.empty:
            continue
        used.add(((dataset, backbone), region))
        for _, r in cell.iterrows():
            y = arms.index(r["arm"])
            color = color_of(r)
            if r["n_pairs"]:
                ax.errorbar(
                    r["mean"],
                    y,
                    xerr=[[r["mean"] - r["ci_lo"]], [r["ci_hi"] - r["mean"]]],
                    fmt="o",
                    color=color,
                    ms=3,
                    capsize=2,
                    lw=1,
                )
            if note := notes(r):
                ax.text(
                    0.98,
                    y,
                    note,
                    transform=ax.get_yaxis_transform(),
                    fontsize=6,
                    ha="right",
                    va="center",
                    color=color,
                )
        ax.axvline(0, color="0.7", lw=0.8)
        for x in rule_lines:
            ax.axvline(x, color="0.4", lw=0.8, ls="--")
        suffix = " (cross-check)" if dataset != "gg2" else ""
        ax.set_title(f"{dataset} {backbone} {region}{suffix}", fontsize=8)
        ax.set_yticks(range(len(arms)), arms, fontsize=7)
        ax.set_xlabel(xlabel, fontsize=7)
        ax.set_xscale("symlog", linthresh=1)
        ax.margins(x=0.08)
        ax.tick_params(axis="x", labelsize=7)
    hide_unused(axes, used)
    return fig


def verdict_notes(r):
    if r["verdict"] == "unplaced":
        return f"{r['n'] - r['n_tool_unplaced']}/{r['n']} placed"
    return "no reference" if r["verdict"] == "no-reference" else ""


def verdict_figure(accuracy):
    fig = interval_figure(
        verdict(accuracy),
        VERDICT_ROWS,
        lambda r: VERDICT_COLORS[r["verdict"]],
        "mean edge-distance difference vs sepp-prod (95% CI)",
        verdict_notes,
        rule_lines=(MARGIN,),
    )
    side_legend(
        fig,
        [
            plt.Line2D([], [], color=c, marker="o", ls="", label=v)
            for v, c in VERDICT_COLORS.items()
        ],
    )
    footnote(
        fig,
        f"Dashed line: pass margin {MARGIN} (upper CI bound at or below it passes). "
        "unplaced: some query is unplaced (k/n placed shown); "
        "no reference: sepp-prod left queries unplaced.",
        "Arms are paired with the sepp-prod cell of the same region and backbone; "
        "oracle excluded.",
        APPLES_CAVEAT,
    )
    return fig


def oracle_gap_figure(accuracy):
    fig = interval_figure(
        oracle_gap(accuracy).assign(verdict="oracle"),
        [("gg2", "B10k")],
        lambda r: "tab:blue",
        "oracle minus hmmalign edge distance (95% CI)",
        lambda r: "",
    )
    footnote(
        fig,
        "Negative: the oracle alignment places closer to truth than hmmalign.",
        APPLES_CAVEAT,
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description="Accuracy figures from accuracy.tsv.")
    parser.add_argument("--accuracy", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    accuracy = pd.read_csv(args.accuracy, sep="\t")
    summary = bin_summary(accuracy)
    both = ["B10k", "Bfull"]
    save(
        metric_figure(summary, "mean_edge_dist", "mean edge distance (placed)", both),
        args.outdir,
        "accuracy_edge_dist.png",
    )
    save(
        metric_figure(summary, "exact_frac", "exact fraction", both),
        args.outdir,
        "accuracy_exact.png",
    )
    save(
        metric_figure(summary, "mean_delta_err", "mean delta error (B10k)", ["B10k"]),
        args.outdir,
        "accuracy_delta_err.png",
    )
    save(verdict_figure(accuracy), args.outdir, "verdict.png")
    save(oracle_gap_figure(accuracy), args.outdir, "oracle_gap.png")


if __name__ == "__main__":
    main()
