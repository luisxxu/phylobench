"""Markdown tables for results.md, built only from the collected result tables."""

import re

import pandas as pd

from phylobench.summaries import (
    BIN_ORDER,
    REFERENCE_ARM,
    bin_summary,
    cell_distribution,
    end_to_end,
    oracle_gap,
    split_run_id,
    verdict,
)

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
]
BACKBONES = ["B10k", "Bfull"]
MIN_REGION_N = 200
SMALL_MARK = "†"
FOCUS_REGIONS = ["full-length", "V4-150"]
LICENSE_COLUMNS = 5
BLOCK = re.compile(r"<!-- table:(?P<key>[\w-]+) -->\n.*?\n<!-- /table -->", re.DOTALL)
CPU_NAME = re.compile(r"(Xeon|EPYC)\D*(\d{4})")


def ordered(values, order):
    return [v for v in order if v in set(values)]


def markdown_table(df: pd.DataFrame, index: bool = True) -> str:
    if index:
        df = df.reset_index()

    def cell(value):
        return "" if pd.isna(value) else str(value).replace("|", "\\|")

    rows = [[cell(c) for c in df.columns], *df.map(cell).values.tolist()]
    lines = [f"| {' | '.join(row)} |" for row in rows]
    lines.insert(1, f"|{'---|' * len(df.columns)}")
    return "\n".join(lines)


def ci_reading(ci_lo: float, ci_hi: float) -> str:
    if ci_lo > 0:
        return "worse"
    if ci_hi < 0:
        return "better"
    return "spans 0"


def effect(row) -> str:
    return f"{row['mean']:+.2f} [{row['ci_lo']:+.2f}, {row['ci_hi']:+.2f}]"


def verdict_cell(row) -> str:
    if pd.isna(row["mean"]):
        return row["verdict"]
    reading = ci_reading(row["ci_lo"], row["ci_hi"])
    return f"{row['verdict']} {effect(row)} {reading}"


def verdict_table(v: pd.DataFrame, dataset: str, backbone: str) -> pd.DataFrame:
    sub = v[(v["dataset"] == dataset) & (v["backbone"] == backbone)]
    regions = ordered(sub["region"], REGION_ORDER)
    cells = sub.assign(cell=sub.apply(verdict_cell, axis=1)).pivot(
        index="arm", columns="region", values="cell"
    )
    table = cells.reindex(index=ordered(sub["arm"], ARM_ORDER), columns=regions)
    table = table.fillna("")
    table.loc["n"] = sub.groupby("region")["n"].first().reindex(regions).astype(str)
    table.index.name = "arm"
    return table


def pass_count_table(v: pd.DataFrame) -> pd.DataFrame:
    rows = {}
    for (dataset, backbone), g in v.groupby(["dataset", "backbone"]):
        for arm, a in g.groupby("arm"):
            large = a[a["n"] >= MIN_REGION_N]
            rows.setdefault(arm, {})[f"{dataset} {backbone}"] = (
                f"{(a['verdict'] == 'pass').sum()}/{len(a)} "
                f"({(large['verdict'] == 'pass').sum()}/{len(large)})"
            )
    table = pd.DataFrame(rows).T.fillna("")
    table = table.reindex(ordered(table.index, ARM_ORDER))
    table.index.name = "arm"
    return table


def fmt(value: float, spec: str = ".2f") -> str:
    return "" if pd.isna(value) else format(value, spec)


def distribution_table(
    dist: pd.DataFrame, dataset: str, regions: list[str]
) -> pd.DataFrame:
    sub = dist[(dist["dataset"] == dataset) & (dist["aln"] != "oracle")]
    rows = []
    for backbone in BACKBONES:
        for region in ordered(regions, REGION_ORDER):
            cell = sub[(sub["backbone"] == backbone) & (sub["region"] == region)]
            for arm in ordered(cell["arm"], ARM_ORDER):
                r = cell[cell["arm"] == arm].iloc[0]
                all_unplaced = r["n_unplaced"] == r["n"]
                rows.append(
                    {
                        "backbone": backbone,
                        "region": region,
                        "arm": arm,
                        "queries": str(r["n"]),
                        "unplaced": str(r["n_unplaced"]),
                        "mean": "unplaced" if all_unplaced else fmt(r["mean"]),
                        "median": "" if all_unplaced else fmt(r["median"], ".1f"),
                        "fraction at 0": "" if all_unplaced else fmt(r["frac_zero"]),
                    }
                )
    return pd.DataFrame(rows)


def bin_table(
    bins: pd.DataFrame, dataset: str, backbone: str, region: str
) -> pd.DataFrame:
    sub = bins[
        (bins["dataset"] == dataset)
        & (bins["backbone"] == backbone)
        & (bins["region"] == region)
        & (bins["aln"] != "oracle")
    ]

    def cell(row):
        if row["n_unplaced"] == row["n"]:
            return "unplaced"
        return f"{row['mean_edge_dist']:.2f}{SMALL_MARK if row['small'] else ''}"

    table = sub.assign(cell=sub.apply(cell, axis=1)).pivot(
        index="arm", columns="bin", values="cell"
    )
    columns = ordered(sub["bin"], BIN_ORDER)
    table = table.reindex(index=ordered(sub["arm"], ARM_ORDER), columns=columns).fillna(
        ""
    )
    table.loc["n"] = (
        sub.groupby("bin", observed=True)["n"].max().reindex(columns).astype(str)
    )
    table.index.name = "arm"
    return table


def oracle_table(gap: pd.DataFrame, dataset: str, backbone: str) -> pd.DataFrame:
    sub = gap[(gap["dataset"] == dataset) & (gap["backbone"] == backbone)]
    regions = ordered(sub["region"], REGION_ORDER)
    cells = sub.assign(cell=sub.apply(effect, axis=1)).pivot(
        index="arm", columns="region", values="cell"
    )
    table = cells.reindex(index=ordered(sub["arm"], ARM_ORDER), columns=regions)
    table.index.name = "arm"
    return table


def short_cpu(cpu_model: str) -> str:
    match = CPU_NAME.search(cpu_model)
    return f"{match[1]} {match[2]}" if match else cpu_model


def e2e_status(row) -> str:
    if row["e2e_status"] == "ok":
        return "ok"
    return "DNF (convert)" if row["convert_failed"] else "DNF"


def cost_table(
    e2e: pd.DataFrame, dataset: str, backbone: str, region: str
) -> pd.DataFrame:
    sub = e2e[
        (e2e["dataset"] == dataset)
        & (e2e["backbone"] == backbone)
        & (e2e["region"] == region)
    ]
    sub = sub.set_index("arm").loc[ordered(sub["arm"], ARM_ORDER)]
    ref = sub.loc[REFERENCE_ARM]
    e2e_ok = sub["e2e_status"] == "ok"
    place_ok = sub["place_status"] == "ok"
    return pd.DataFrame(
        {
            "aln": sub["aln"],
            "queries": sub["query_n"].map(lambda n: f"{n:.0f}"),
            "place s": sub["place_s"].map(lambda s: fmt(s, ".1f")),
            "e2e s": sub["e2e_s"].map(lambda s: fmt(s, ".1f")),
            "e2e vs sepp-prod": (sub["e2e_s"] / ref["e2e_s"]).where(e2e_ok).map(fmt),
            "peak mem MB": sub["peak_mem_mb"].map(lambda m: fmt(m, ".0f")),
            "mem vs sepp-prod": (sub["peak_mem_mb"] / ref["peak_mem_mb"])
            .where(place_ok)
            .map(fmt),
            "gpu peak MB": sub["gpu_peak_mb"].map(lambda m: fmt(m, ".0f")),
            "mem source": sub["mem_source"],
            "cpu": sub["cpu_model"].map(short_cpu),
            "e2e status": sub.apply(e2e_status, axis=1),
        }
    ).reset_index()


def runs_overview(runs: pd.DataFrame) -> pd.DataFrame:
    grouped = runs.groupby("phase")
    counts = {
        "rows": grouped.size(),
        "ok": grouped["status"].apply(lambda s: (s == "ok").sum()),
        "DNF": grouped["status"].apply(lambda s: (s == "DNF").sum()),
        "cgroup": grouped["mem_source"].apply(lambda s: (s == "cgroup").sum()),
        "psutil": grouped["mem_source"].apply(lambda s: (s == "psutil").sum()),
    }
    return pd.DataFrame({k: v.astype(str) for k, v in counts.items()}).reset_index()


def hardware_table(e2e: pd.DataFrame) -> pd.DataFrame:
    machine = e2e["backbone"] + " " + e2e["cpu_model"].map(short_cpu)
    counts = e2e.groupby(["arm", machine]).size().unstack(fill_value=0)
    counts = counts.reindex(ordered(counts.index, ARM_ORDER))
    counts.columns.name = None
    counts.index.name = "arm"
    return counts.astype(str)


def concordance_table(concordance: pd.DataFrame) -> pd.DataFrame:
    df = split_run_id(concordance)
    df["at_zero"] = df["edge_dist"].eq(0)
    grouped = df.groupby(["dataset", "backbone", "arm"])
    n = grouped.size()
    unplaced = n - grouped["edge_dist"].count()
    out = pd.DataFrame(
        {
            "queries": n.astype(str),
            "unplaced": unplaced.astype(int).astype(str),
            "mean": grouped["edge_dist"].mean().map(fmt),
            "median": grouped["edge_dist"].median().map(lambda m: fmt(m, ".1f")),
            "fraction at 0": (grouped["at_zero"].sum() / (n - unplaced)).map(fmt),
        }
    ).reset_index()
    out["arm_rank"] = out["arm"].map(ARM_ORDER.index)
    out["backbone_rank"] = out["backbone"].map(BACKBONES.index)
    out = out.sort_values(["dataset", "backbone_rank", "arm_rank"])
    return out.drop(columns=["arm_rank", "backbone_rank"]).reset_index(drop=True)


def license_table(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("| Arm / tool"))
    rows = []
    for line in lines[start:]:
        if not line.startswith("|"):
            break
        rows.append([c.strip() for c in line.strip("|").split(" | ")][:LICENSE_COLUMNS])
    return pd.DataFrame(rows[2:], columns=rows[0])


VERDICT_NOTE = (
    "Cell: pre-registered verdict (pass if the upper 95% bootstrap CI of the mean "
    "paired difference vs sepp-prod is at most 0.25 edges), mean difference "
    "[95% CI] in edges (tool minus sepp-prod; negative = tool closer), then whether "
    "the CI is above 0 (`worse`), below 0 (`better`) or spans 0. `n` = queries scored."
)
BIN_NOTE = (
    "Mean edge distance per similarity bin; `†` = fewer than 100 queries; `n` = "
    "queries in the bin. Descriptive only, never pooled across backbones."
)
COST_NOTE = (
    "`place s` = placement wall time; `e2e s` = align + place + convert; ratios are "
    "to sepp-prod in the same cell and blank where the cell did not finish. Peak "
    "memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only)."
)
BACKBONE_CELLS = [("gg2", "B10k"), ("gg2", "Bfull"), ("sim10k", "B10k")]


def build_blocks(
    runs: pd.DataFrame,
    accuracy: pd.DataFrame,
    concordance: pd.DataFrame | None,
    licenses_text: str,
) -> dict[str, str]:
    """Table key -> markdown block; numbered by position, empty tables skipped."""
    v = verdict(accuracy)
    bins = bin_summary(accuracy)
    e2e = end_to_end(runs)
    specs = [
        ("runs-overview", "Rows of runs.tsv per phase.", runs_overview(runs), False),
        (
            "hardware",
            "Placed cells per arm, backbone and CPU.",
            hardware_table(e2e),
            True,
        ),
        (
            "pass-counts",
            "Regions passed per arm: all regions (regions with n >= 200).",
            pass_count_table(v),
            True,
        ),
    ]
    for dataset, backbone in BACKBONE_CELLS:
        caption = f"Verdict vs sepp-prod, {dataset} at {backbone}\\*. {VERDICT_NOTE}"
        table = verdict_table(v, dataset, backbone)
        specs.append((f"verdict-{dataset}-{backbone}", caption, table, True))
    specs.append(
        (
            "distribution-gg2",
            "Edge distance per cell, gg2 full-length and V4-150 (sepp-prod included): "
            "unplaced queries, mean and median over placed queries, fraction of placed "
            "queries at distance 0.",
            distribution_table(cell_distribution(accuracy), "gg2", FOCUS_REGIONS),
            False,
        )
    )
    specs.append(
        (
            "oracle-gg2-B10k",
            "Oracle minus hmmalign alignment, same arm, gg2 at B10k\\*: mean "
            "edge-distance difference [95% CI]; negative = oracle alignment is closer.",
            oracle_table(oracle_gap(accuracy), "gg2", "B10k"),
            True,
        )
    )
    for region in FOCUS_REGIONS:
        for backbone in BACKBONES:
            caption = (
                f"Per-bin mean edge distance, gg2 {region} at {backbone}\\*. {BIN_NOTE}"
            )
            table = bin_table(bins, "gg2", backbone, region)
            specs.append((f"bins-gg2-{backbone}-{region}", caption, table, True))
    cost_cells = [("gg2", region) for region in FOCUS_REGIONS]
    cost_cells += [("emp-10k", "V4-150"), ("emp-100k", "V4-150")]
    for dataset, region in cost_cells:
        for backbone in BACKBONES:
            caption = f"Cost, {dataset} {region} at {backbone}\\*. {COST_NOTE}"
            table = cost_table(e2e, dataset, backbone, region)
            key = f"cost-{dataset}-{backbone}-{region}"
            specs.append((key, caption, table, False))
    if concordance is not None:
        specs.append(
            (
                "concordance",
                "EMP concordance with sepp-prod (edge distance between each arm's "
                "placement and sepp-prod's; not accuracy). Mean and median over "
                "placed ASVs.",
                concordance_table(concordance),
                False,
            )
        )
    specs.append(
        (
            "licenses",
            "Licenses (flag only; from LICENSES.md).",
            license_table(licenses_text),
            False,
        )
    )
    return {
        key: (
            f"<!-- table:{key} -->\n**Table {i}.** {caption}\n\n"
            f"{markdown_table(table, index)}\n<!-- /table -->"
        )
        for i, (key, caption, table, index) in enumerate(specs, 1)
        if len(table)
    }


def extract_blocks(text: str) -> dict[str, str]:
    return {m["key"]: m[0] for m in BLOCK.finditer(text)}


def paste_blocks(text: str, blocks: dict[str, str]) -> str:
    for key, block in blocks.items():
        if f"<!-- table:{key} -->" not in text:
            raise ValueError(f"results.md has no <!-- table:{key} --> marker")
        text = BLOCK.sub(lambda m: block if m["key"] == key else m[0], text)
    return text
