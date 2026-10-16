"""Tests for RAxML-8 info to EPA-ng model conversion.

The GG2 fixture tests/fixtures/gg2-2022.10/raxml-info.txt is the byte-identical
member 9e8f6bb7-1444-4e32-b5fd-92cd5823d3ee/data/raxml-info.txt
(sha256 1ad6d7b1b2c696a93ff5dcf2be1e677817bb50063474e41c0c0d3dfc22449f76) of the
GG2 2022.10 SEPP reference artifact sepp-ref.qza
(sha256 fa5d32f936a39735a4f99552044ed4fcbae0d91fb48a679f9d329e2a708eb973).
"""

import subprocess
import sys
from pathlib import Path

import pytest

from phylobench.models import GtrModel, epa_ng_model, parse_raxml_info

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"
SCRIPT = ROOT / "workflow" / "scripts" / "write_epa_ng_model.py"
TINY_INFO = FIXTURES / "tiny" / "raxml-info.txt"
GG2_INFO = FIXTURES / "gg2-2022.10" / "raxml-info.txt"

LABELS = "ac ag at cg ct gt"
RATES = "1.0 2.0 3.0 4.0 5.0 6.0"
FREQS = "0.1 0.2 0.3 0.4"


def info(rates=RATES, freqs=FREQS, labels=LABELS, extra=""):
    return (
        f"Base frequencies: {freqs}\n"
        f"alpha[0]: 1.000000 rates[0] {labels}: {rates}\n{extra}"
    )


def test_gg2_info_to_epa_ng_model():
    assert epa_ng_model(parse_raxml_info(GG2_INFO.read_text())) == (
        "GTR{0.900673/2.417312/1.491932/0.621988/3.987304/1.0}"
        "+FU{0.262062/0.222878/0.306557/0.208502}+G4m{1.0}"
    )


def test_tiny_info_to_epa_ng_model():
    assert epa_ng_model(parse_raxml_info(TINY_INFO.read_text())) == (
        "GTR{0.947786/2.643649/1.480482/0.78679/4.471837/1.0}"
        "+FU{0.250954/0.234726/0.31899/0.19533}+G4m{1.0}"
    )


def test_rates_map_positionally_to_ac_ag_at_cg_ct_gt():
    model = parse_raxml_info(info())
    assert model == GtrModel((1.0, 2.0, 3.0, 4.0, 5.0, 6.0), (0.1, 0.2, 0.3, 0.4), 1.0)
    assert epa_ng_model(model).startswith("GTR{1.0/2.0/3.0/4.0/5.0/6.0}")


def test_permuted_rate_labels_rejected():
    with pytest.raises(ValueError, match="rate labels"):
        parse_raxml_info(info(labels="ac ag cg at ct gt"))


@pytest.mark.parametrize("bad", ["nan", "inf", "-1", "0"])
def test_rejects_bad_numbers(bad):
    with pytest.raises(ValueError):
        parse_raxml_info(info(rates=f"{bad} 2 3 4 5 6"))
    with pytest.raises(ValueError):
        parse_raxml_info(info(freqs=f"{bad} 0.2 0.3 0.4"))


def test_rejects_wrong_counts():
    with pytest.raises(ValueError):
        parse_raxml_info(info(rates="1 2 3 4 5"))
    with pytest.raises(ValueError):
        parse_raxml_info(info(freqs="0.3 0.3 0.4"))


def test_rejects_multiple_partitions():
    extra = f"alpha[1]: 1.000000 rates[1] {LABELS}: {RATES}\n"
    with pytest.raises(ValueError):
        parse_raxml_info(info(extra=extra))


def test_rejects_missing_lines():
    with pytest.raises(ValueError):
        parse_raxml_info(f"alpha[0]: 1.0 rates[0] {LABELS}: {RATES}\n")
    with pytest.raises(ValueError):
        parse_raxml_info(f"Base frequencies: {FREQS}\n")


def test_write_epa_ng_model_script(tmp_path):
    aln = tmp_path / "aln.fasta"
    aln.write_text(">a\nACGT\nAC\n>b\nACGTACGTAC\n")
    out = tmp_path / "model"
    subprocess.run(
        [sys.executable, SCRIPT, TINY_INFO, aln, out],
        check=True,
        env={"PYTHONPATH": str(ROOT / "src")},
    )
    assert out.read_text() == (
        "GTR{0.947786/2.643649/1.480482/0.78679/4.471837/1.0}"
        "+FU{0.250954/0.234726/0.31899/0.19533}+G4m{1.0}, noname = 1-6\n"
    )
