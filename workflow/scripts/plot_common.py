import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REFERENCE_ARM = "sepp-prod"
REGION_ORDER = [
    "full-length",
    "V4-150",
    "V3-V4",
    "V1-V2",
    "V1-V3",
    "V4-V5",
    "V5-V6",
    "V6-V8",
    "V7-V9",
    *(
        f"{r}c"
        for r in (
            "V4-150",
            "V3-V4",
            "V1-V2",
            "V1-V3",
            "V4-V5",
            "V5-V6",
            "V6-V8",
            "V7-V9",
        )
    ),
]
ARM_ORDER = [
    REFERENCE_ARM,
    "sepp-pinned",
    "pplacer",
    "pplacer-scampp",
    "epa-ng",
    "epa-ng-scampp",
    "bscampp",
    "apples2",
    "dipper-m",
    "dipper-m-gpu",
    "dipper-r",
]
DPI = 200
APPLES_CAVEAT = (
    "GG2 truth was seeded by uDance (APPLES-2): apples2 accuracy is optimistic."
)


def arm_color(arm: str):
    if arm == REFERENCE_ARM:
        return "black"
    others = [a for a in ARM_ORDER if a != REFERENCE_ARM]
    return plt.get_cmap("tab10")(others.index(arm) % 10)


def present_arms(arms) -> list[str]:
    return [a for a in ARM_ORDER if a in set(arms)]


def present_regions(regions) -> list[str]:
    return [r for r in REGION_ORDER if r in set(regions)]


def short_cpu(cpu_model: str) -> str:
    match = re.search(r"(Xeon|EPYC)\D*(\d{4})", cpu_model)
    return f"{match[1]} {match[2]}" if match else cpu_model


def arm_labels(runs: pd.DataFrame) -> dict[str, str]:
    """Arm names, tagged with the CPU when an arm ran on other than the majority CPU."""
    timed = runs[runs["arm"].notna() & (runs["phase"] != "evaluate")]
    majority = timed["cpu_model"].mode().iat[0]
    off = timed[timed["cpu_model"] != majority].groupby("arm")["cpu_model"].first()
    return {
        arm: f"{arm} ({short_cpu(off[arm])})" if arm in off.index else arm
        for arm in timed["arm"].unique()
    }


def hardware_note(runs: pd.DataFrame) -> str | None:
    labels = arm_labels(runs)
    tagged = sorted(a for a, label in labels.items() if label != a)
    if not tagged:
        return None
    return (
        f"Hardware: {', '.join(tagged)} ran on a different CPU than the other arms "
        "(tag in legend); runtimes across hardware are not strictly comparable."
    )


def grid(row_keys, col_keys, *, panel=(3.0, 2.6)):
    fig, axes = plt.subplots(
        len(row_keys),
        len(col_keys),
        figsize=(panel[0] * len(col_keys) + 2.6, panel[1] * len(row_keys) + 1.2),
        squeeze=False,
        sharey="row",
        layout="constrained",
    )
    return fig, {
        (r, c): axes[i][j]
        for i, r in enumerate(row_keys)
        for j, c in enumerate(col_keys)
    }


def region_grid(regions_by_key, *, per_row=5, panel=(3.0, 2.6)):
    chunks = {
        (key, i // per_row): regions[i : i + per_row]
        for key, regions in regions_by_key.items()
        for i in range(0, len(regions), per_row)
    }
    fig, panels = grid(list(chunks), range(per_row), panel=panel)
    axes = {}
    for ((key, i), j), ax in panels.items():
        if j < len(chunks[(key, i)]):
            axes[(key, chunks[(key, i)][j])] = ax
        else:
            ax.set_visible(False)
    return fig, axes


def hide_unused(axes, used):
    for key, ax in axes.items():
        if key not in used:
            ax.set_visible(False)


def side_legend(fig, handles):
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=8,
        frameon=False,
    )


def legend(fig, arms, labels):
    handles = [
        plt.Line2D(
            [], [], color=arm_color(a), marker="o", label=labels.get(a, a), lw=1.5
        )
        for a in arms
    ]
    side_legend(fig, handles)


def footnote(fig, *lines):
    width = int(fig.get_figwidth() * 17)
    wrapped = [textwrap.fill(line, width) for line in lines if line]
    fig.text(0.0, -0.01, "\n".join(wrapped), fontsize=7, va="top")


def save(fig, outdir: Path, name: str):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
