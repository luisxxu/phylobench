"""Primer-matched amplicon regions and helpers for the shared hmmalign alignment."""

import re
from collections import Counter

GAPS = "-."
IUPAC = {
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}  # fmt: skip
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
PRIMERS = {
    "V4-150": (
        "GTGYCAGCMGCCGCGGTAA",
        "GGACTACNVGGGTWTCTAAT",
    ),  # 515F Parada, 806R Apprill
    "V3-V4": ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"),  # 341F, 805R
    "V1-V2": ("AGAGTTTGATCMTGGCTCAG", "GCTGCCTCCCGTAGGAGT"),  # 27F, 338R
    "V5-V6": ("AGGATTAGATACCCTGGTA", "CGACRRCCATGCANCACCT"),  # 784F, 1064R
}
# V4-150 is the first 150 nt after 515F, but 806R must still be present.
TRIM = {"V4-150": 150}
REGIONS = ("full-length", *PRIMERS)

_DELETE_INSERTS = str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz.")
_DELETE_GAPS = str.maketrans("", "", GAPS)


def read_fasta(path):
    """Yield (name, sequence) pairs; multi-line records are joined."""
    name, parts = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def primer_pattern(primer):
    return re.compile("".join(IUPAC.get(c, c) for c in primer))


def ungap(seq):
    return seq.translate(_DELETE_GAPS)


def amplicon_span(seq, region):
    """(start, end) on the ungapped sequence, primers excluded, or a reason string."""
    seq = seq.upper()
    if region == "full-length":
        return 0, len(seq)
    forward, reverse = PRIMERS[region]
    forward_hit = primer_pattern(forward).search(seq)
    if forward_hit is None:
        return "no forward primer"
    start = forward_hit.end()
    reverse_hit = primer_pattern(reverse_complement(reverse)).search(seq, start)
    if reverse_hit is None:
        return "no reverse primer"
    end = start + TRIM[region] if region in TRIM else reverse_hit.start()
    if end > reverse_hit.start() or end <= start:
        return "too short"
    return start, end


def split_spans(pairs, region):
    """Return ({name: span}, {name: exclusion reason}) for (name, sequence) pairs."""
    spans, exclusions = {}, {}
    for name, seq in pairs:
        result = amplicon_span(seq, region)
        if isinstance(result, str):
            exclusions[name] = result
        else:
            spans[name] = result
    return spans, exclusions


def exclusion_counts(exclusions):
    return dict(Counter(exclusions.values()))


def project(span, aligned_row):
    """Map an ungapped half-open span onto half-open columns of the aligned row."""
    columns = [i for i, c in enumerate(aligned_row) if c not in GAPS]
    start, end = span
    return columns[start], columns[end - 1] + 1


def check_ungaps_to(aligned, full_length):
    if set(aligned) != set(full_length):
        raise ValueError("aligned and full-length FASTA have different names")
    for name, row in aligned.items():
        if ungap(row).upper() != full_length[name].upper():
            raise ValueError(
                f"{name}: aligned row does not ungap to full-length sequence"
            )


def stockholm_lines(rows):
    """One-block Stockholm with an all-x RF line (hmmbuild --hand keeps all columns)."""
    width = None
    yield "# STOCKHOLM 1.0"
    for name, row in rows:
        width = len(row) if width is None else width
        if len(row) != width:
            raise ValueError(f"{name}: row width {len(row)} differs from {width}")
        yield f"{name} {row}"
    yield f"#=GC RF {'x' * (width or 0)}"
    yield "//"


def mask_inserts(row):
    """Drop hmmalign insert columns (lowercase residues and '.' gaps)."""
    return row.translate(_DELETE_INSERTS)
