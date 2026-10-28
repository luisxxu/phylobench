"""Per-bin accuracy summaries, paired differences vs sepp-prod, and end-to-end cost."""

import numpy as np
import pandas as pd

CELL_KEYS = ["dataset", "backbone", "region", "aln", "arm"]
BIN_ORDER = ["<75", "75-80", "80-85", "85-90", "90-95", "95-97", "97-99", "99-99.9"]
BIN_ORDER.append(">=99.9")
REFERENCE_ARM = "sepp-prod"
MARGIN = 0.25
SMALL_BIN_N = 100
N_BOOT = 10_000
PLACE_PHASES = ("place", "align+place", "align+place+convert")
BOOT_CHUNK = 1_000
UNPAIRED = {"n_pairs": 0, "n_tool_unplaced": 0, "n_ref_unplaced": 0}
NO_STATS = {"mean": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}


def split_run_id(df: pd.DataFrame) -> pd.DataFrame:
    parts = df["run_id"].str.split("/", expand=True)
    if parts.shape[1] != len(CELL_KEYS) or parts.isna().any(axis=None):
        raise ValueError("run_id must be dataset/backbone/region/aln/arm")
    parts.columns = CELL_KEYS
    return df.join(parts)


def bin_summary(accuracy: pd.DataFrame) -> pd.DataFrame:
    df = split_run_id(accuracy)
    df["bin"] = pd.Categorical(df["bin"], BIN_ORDER, ordered=True)
    grouped = df.groupby([*CELL_KEYS, "bin"], observed=True)
    out = pd.DataFrame(
        {
            "n": grouped.size(),
            "n_unplaced": grouped["edge_dist"].apply(lambda s: s.isna().sum()),
            "mean_edge_dist": grouped["edge_dist"].mean(),
            "exact_frac": grouped["exact"].sum() / grouped.size(),
            "mean_delta_err": grouped["delta_err"].mean(),
        }
    )
    out["small"] = out["n"] < SMALL_BIN_N
    return out.reset_index()


def paired_difference(
    tool: pd.DataFrame, ref: pd.DataFrame, *, n_boot: int = N_BOOT, seed: int = 0
) -> dict:
    """tool minus ref edge_dist per query, with a bootstrap 95% CI of the mean.

    Queries unplaced on either side are left out of the mean and counted.
    """
    if set(tool["query"]) != set(ref["query"]):
        raise ValueError("tool and reference score different query sets")
    paired = tool.merge(ref, on="query", suffixes=("_tool", "_ref"))
    tool_nan = paired["edge_dist_tool"].isna()
    ref_nan = paired["edge_dist_ref"].isna()
    diffs = (paired["edge_dist_tool"] - paired["edge_dist_ref"]).dropna().to_numpy()
    counts = {
        "n_pairs": len(diffs),
        "n_tool_unplaced": int((tool_nan & ~ref_nan).sum()),
        "n_ref_unplaced": int(ref_nan.sum()),
    }
    if not len(diffs):
        return {**counts, **NO_STATS}
    rng = np.random.default_rng(seed)
    means = np.concatenate(
        [
            diffs[
                rng.integers(0, len(diffs), (min(BOOT_CHUNK, n_boot - i), len(diffs)))
            ].mean(axis=1)
            for i in range(0, n_boot, BOOT_CHUNK)
        ]
    )
    lo, hi = np.percentile(means, [2.5, 97.5])
    return {**counts, "mean": diffs.mean(), "ci_lo": lo, "ci_hi": hi}


def verdict_state(stats: dict | None, margin: float) -> str:
    if stats is None or stats["n_ref_unplaced"]:
        return "no-reference"
    if stats["n_tool_unplaced"]:
        return "unplaced"
    return "pass" if stats["ci_hi"] <= margin else "fail"


def scored_cells(accuracy: pd.DataFrame) -> dict[tuple, pd.DataFrame]:
    return dict(iter(split_run_id(accuracy).groupby(CELL_KEYS)))


def verdict(accuracy: pd.DataFrame, margin: float = MARGIN) -> pd.DataFrame:
    """Each hmmalign or raw arm against sepp-prod per region."""
    cells = scored_cells(accuracy)
    rows = []
    for (dataset, backbone, region, aln, arm), tool in cells.items():
        if arm == REFERENCE_ARM or aln == "oracle":
            continue
        ref = cells.get((dataset, backbone, region, "raw", REFERENCE_ARM))
        stats = None if ref is None else paired_difference(tool, ref)
        rows.append(
            {
                "dataset": dataset,
                "backbone": backbone,
                "region": region,
                "aln": aln,
                "arm": arm,
                "n": len(tool),
                **(stats or {**UNPAIRED, **NO_STATS}),
                "verdict": verdict_state(stats, margin),
            }
        )
    return pd.DataFrame(rows)


def oracle_gap(accuracy: pd.DataFrame) -> pd.DataFrame:
    """Oracle minus hmmalign edge_dist, same arm; negative means oracle is closer."""
    cells = scored_cells(accuracy)
    rows = []
    for (dataset, backbone, region, aln, arm), oracle in cells.items():
        hmmalign = cells.get((dataset, backbone, region, "hmmalign", arm))
        if aln == "oracle" and hmmalign is not None:
            keys = {"dataset": dataset, "backbone": backbone, "region": region}
            rows.append({**keys, "arm": arm, **paired_difference(oracle, hmmalign)})
    return pd.DataFrame(rows)


def end_to_end(runs: pd.DataFrame) -> pd.DataFrame:
    """One row per placed cell: place-only wall vs the sum of align + place + convert.

    ref-prep and per-backbone hmmbuild are one-time costs and stay out. A failed convert
    marks e2e DNF while the placement and its accuracy rows stay valid.
    """
    cells = runs[runs["arm"].notna() & (runs["aln"] != "oracle")]
    out = cells[cells["phase"].isin(PLACE_PHASES)].rename(
        columns={"status": "place_status"}
    )
    align_keys = ["dataset", "backbone", "region", "aln"]
    align = runs[
        runs["arm"].isna() & runs["region"].notna() & (runs["phase"] == "align")
    ]
    align = align[[*align_keys, "wall_s", "status"]].rename(
        columns={"wall_s": "align_s", "status": "align_status"}
    )
    convert = cells[cells["phase"] == "convert"][[*CELL_KEYS, "wall_s", "status"]]
    convert = convert.rename(
        columns={"wall_s": "convert_s", "status": "convert_status"}
    )

    out = out.merge(align, how="left", on=align_keys).merge(
        convert, how="left", on=CELL_KEYS
    )
    needs_align = (out["phase"] == "place") & (out["aln"] == "hmmalign")
    out["align_s"] = out["align_s"].where(needs_align, 0.0)
    out["align_status"] = out["align_status"].where(needs_align, "ok")
    out["convert_s"] = out["convert_s"].fillna(0.0)
    out["convert_failed"] = out["convert_status"] == "DNF"
    out["place_s"] = out["wall_s"]
    out["e2e_s"] = out["align_s"] + out["wall_s"] + out["convert_s"]
    parts_ok = (
        (out["place_status"] == "ok")
        & (out["align_status"] == "ok")
        & ~out["convert_failed"]
    )
    out["e2e_status"] = np.where(parts_ok, "ok", "DNF")
    columns = [
        *CELL_KEYS,
        "query_n",
        "place_s",
        "e2e_s",
        "peak_mem_mb",
        "mem_source",
        "gpu_peak_mb",
        "cpu_model",
        "place_status",
        "convert_failed",
        "e2e_status",
    ]
    return out[columns].reset_index(drop=True)
