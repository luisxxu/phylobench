import math
import re
from typing import NamedTuple

RATE_LABELS = ("ac", "ag", "at", "cg", "ct", "gt")

_FREQS_LINE = re.compile(r"^Base frequencies: (.+)$", re.M)
_ALPHA_LINE = re.compile(r"^alpha\[0\]: (\S+) rates\[0\] (.+?): (.+)$", re.M)
_EXTRA_PARTITION = re.compile(r"^alpha\[[1-9]\d*\]", re.M)


class GtrModel(NamedTuple):
    rates: tuple[float, ...]
    freqs: tuple[float, ...]
    alpha: float


def _numbers(text: str, count: int, what: str) -> tuple[float, ...]:
    values = tuple(float(x) for x in text.split())
    if len(values) != count:
        raise ValueError(f"expected {count} {what}, got {len(values)}")
    if not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError(f"{what} must be finite and positive: {values}")
    return values


def parse_raxml_info(text: str) -> GtrModel:
    freq_lines = _FREQS_LINE.findall(text)
    alpha_lines = _ALPHA_LINE.findall(text)
    if len(freq_lines) != 1 or len(alpha_lines) != 1:
        raise ValueError(
            "expected exactly one 'Base frequencies' and one 'alpha[0]' line"
        )
    if _EXTRA_PARTITION.search(text):
        raise ValueError("multiple partitions are not supported")
    alpha, labels, rates = alpha_lines[0]
    if tuple(labels.split()) != RATE_LABELS:
        raise ValueError(f"rate labels must be {' '.join(RATE_LABELS)}, got {labels!r}")
    return GtrModel(
        rates=_numbers(rates, 6, "rates"),
        freqs=_numbers(freq_lines[0], 4, "base frequencies"),
        alpha=_numbers(alpha, 1, "alpha")[0],
    )


def epa_ng_model(model: GtrModel) -> str:
    rates = "/".join(map(repr, model.rates))
    freqs = "/".join(map(repr, model.freqs))
    return f"GTR{{{rates}}}+FU{{{freqs}}}+G4m{{{model.alpha!r}}}"
