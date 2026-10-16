import time

import numpy as np

from phylobench.similarity import similarity

TIME_LIMIT_S = 300


def random_codes(prefix: str, n: int, n_columns: int, rng) -> dict[str, str]:
    alphabet = np.frombuffer(b"ACGT", dtype=np.uint8)
    rows = alphabet[rng.integers(0, 4, size=(n, n_columns))]
    return {f"{prefix}{i}": row.tobytes().decode() for i, row in enumerate(rows)}


def test_10k_by_10k_by_1500_columns():
    rng = np.random.default_rng(0)
    refs = random_codes("r", 10_000, 1500, rng)
    queries = random_codes("q", 10_000, 1500, rng)

    start = time.perf_counter()
    df = similarity(queries, refs)
    elapsed = time.perf_counter() - start

    print(f"10k x 10k x 1500 similarity: {elapsed:.1f} s")
    assert len(df) == 10_000
    assert elapsed < TIME_LIMIT_S
