import sys
from pathlib import Path

from phylobench.models import epa_ng_model, parse_raxml_info


def first_record_width(fasta: Path) -> int:
    width = 0
    with fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if width:
                    break
            else:
                width += len(line.strip())
    return width


info, aln, out = map(Path, sys.argv[1:4])
model = epa_ng_model(parse_raxml_info(info.read_text()))
out.write_text(f"{model}, noname = 1-{first_record_width(aln)}\n")
