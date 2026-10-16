"""Query-to-backbone similarity and similarity bins.

sim(q) = max over reference leaves of matches / n_columns, where n_columns is the
query's non-gap columns inside its own window. Only unambiguous A/C/G/T (U read as T,
any case) can match; a reference gap or ambiguity code is a mismatch, and an ambiguous
query base still counts in n_columns. Memory is about 16 * R * L bytes for R reference
leaves and L alignment columns (float32 one-hot of the references).
"""

import numpy as np
import pandas as pd

BIN_FLOORS_PER_MILLE = (750, 800, 850, 900, 950, 970, 990, 999)
BIN_LABELS = ("<75", "75-80", "80-85", "85-90", "90-95", "95-97", "97-99", "99-99.9")
TOP_BIN = ">=99.9"
MIN_BIN_SIZE = 100

_GAP = 4
_AMBIGUOUS = 5


def _code_table() -> np.ndarray:
    table = np.full(256, _AMBIGUOUS, dtype=np.uint8)
    for base, code in zip("ACGT", range(4), strict=True):
        table[ord(base)] = table[ord(base.lower())] = code
    table[ord("U")] = table[ord("u")] = 3
    table[ord("-")] = table[ord(".")] = _GAP
    return table


_CODES = _code_table()


def _encode(seqs: dict[str, str], n_columns: int) -> np.ndarray:
    codes = np.empty((len(seqs), n_columns), dtype=np.uint8)
    for i, (name, seq) in enumerate(seqs.items()):
        if len(seq) != n_columns:
            raise ValueError(f"{name}: length {len(seq)} != {n_columns} columns")
        codes[i] = _CODES[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    return codes


def _one_hot(codes: np.ndarray) -> np.ndarray:
    out = np.zeros((codes.shape[0], codes.shape[1], 4), dtype=np.float32)
    for base in range(4):
        out[:, :, base] = codes == base
    return out.reshape(codes.shape[0], -1)


def similarity(
    queries: dict[str, str],
    refs: dict[str, str],
    windows: dict[str, tuple[int, int]] | None = None,
    chunk: int = 1000,
) -> pd.DataFrame:
    """One row per query: query, matches, n_columns, sim.

    windows maps a query to half-open column bounds (default: all columns). Float32
    products of 0/1 values sum exactly while the alignment has under 2**24 columns.
    """
    if not refs:
        raise ValueError("no reference sequences")
    n_columns = len(next(iter(refs.values())))
    ref_codes = _encode(refs, n_columns)
    query_codes = _encode(queries, n_columns)
    ref_hot = _one_hot(ref_codes)
    windows = windows or {}

    inside = np.zeros(query_codes.shape, dtype=bool)
    for i, name in enumerate(queries):
        start, end = windows.get(name, (0, n_columns))
        if not 0 <= start <= end <= n_columns:
            raise ValueError(f"{name}: window {(start, end)} outside 0..{n_columns}")
        inside[i, start:end] = True
    counted = inside & (query_codes != _GAP)
    n_counted = counted.sum(axis=1)
    empty = np.flatnonzero(n_counted == 0)
    if empty.size:
        raise ValueError(f"{list(queries)[empty[0]]}: no non-gap columns in its window")

    matches = np.empty(len(queries), dtype=np.int64)
    for lo in range(0, len(queries), chunk):
        hi = lo + chunk
        hot = _one_hot(query_codes[lo:hi]) * np.repeat(inside[lo:hi], 4, axis=1)
        matches[lo:hi] = (hot @ ref_hot.T).max(axis=1)
    return pd.DataFrame(
        {
            "query": list(queries),
            "matches": matches,
            "n_columns": n_counted.astype(np.int64),
            "sim": matches / n_counted,
        }
    )


def similarity_bin(matches: np.ndarray, n_columns: np.ndarray) -> np.ndarray:
    """Bin label per query using integer arithmetic only (no float boundary errors)."""
    matches = np.asarray(matches, dtype=np.int64)
    n_columns = np.asarray(n_columns, dtype=np.int64)
    labels = np.full(matches.shape, BIN_LABELS[0], dtype=object)
    for floor, label in zip(
        BIN_FLOORS_PER_MILLE, BIN_LABELS[1:] + (TOP_BIN,), strict=True
    ):
        labels[matches * 1000 >= floor * n_columns] = label
    return labels.astype(str)


def bin_report(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-bin lineage composition and query-length summary for ONE backbone.

    table has columns query, bin, lineage, length. Never pass rows from more than one
    backbone size: bins are not pooled across backbones.
    """
    composition = table.groupby(["bin", "lineage"]).size().reset_index(name="n")
    summary = (
        table.groupby("bin")["length"]
        .agg(
            n="size",
            length_min="min",
            length_median="median",
            length_max="max",
        )
        .reset_index()
    )
    summary.insert(2, "flagged", summary["n"] < MIN_BIN_SIZE)
    return composition, summary
