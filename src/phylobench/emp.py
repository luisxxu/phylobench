"""EMP Deblur V4 ASV query sets (docs/datasets/emp.md)."""

import re
from collections import Counter

import numpy as np

EMP_QUERY_N = {"emp-10k": 10_000, "emp-100k": 100_000}
EMP_REGION = "V4-150"
EMP_LENGTH = 150
REFERENCE_ARM = "sepp-prod"

_ASV = re.compile(rf"[ACGT]{{{EMP_LENGTH}}}")


def asv_failures(ids: list[str]) -> list[str]:
    """IDs that are not uppercase ACGT of EMP_LENGTH, or that occur more than once."""
    counts = Counter(ids)
    return [i for i in counts if not _ASV.fullmatch(i) or counts[i] > 1]


def query_set(ids: list[str], n: int, seed: int = 0) -> list[str]:
    """First n of a seeded permutation of the sorted IDs, so smaller sets nest."""
    if n > len(ids):
        raise ValueError(f"asked for {n} ASVs but only {len(ids)} are available")
    ordered = sorted(ids)
    order = np.random.default_rng(seed).permutation(len(ordered))
    return [ordered[i] for i in order[:n].tolist()]


def asv_names(n: int) -> list[str]:
    """Short query names: 150-character names are fragile across tools."""
    return [f"asv{i:06d}" for i in range(1, n + 1)]
