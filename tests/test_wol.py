import subprocess
import sys
from pathlib import Path

import treeswift
from test_column_regions import aligned_row
from test_gg2 import lineages_by_cut

from phylobench.regions import read_fasta
from phylobench.wol import (
    clade_capacity,
    dry_run_capacity,
    genome_leaves,
    label_counts,
    length_quantiles,
    lineage_depths,
    neighbour_leaves,
    region_report,
    verdicts,
)

ROOT = Path(__file__).parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "tiny"
SCRIPT = ROOT / "workflow" / "scripts" / "probe_wol.py"


def wol_fixture():
    """The tiny backbone with every other leaf renamed to a G label."""
    tree = treeswift.read_tree_newick((FIXTURE / "tree.nwk").read_text())
    seqs = dict(read_fasta(FIXTURE / "aligned-dna-sequences.fasta"))
    names = {
        leaf.label: f"G{i:09d}" if i % 2 else f"RS_GCF_{i:09d}.1"
        for i, leaf in enumerate(tree.traverse_leaves())
    }
    for leaf in tree.traverse_leaves():
        leaf.label = names[leaf.label]
    return tree, {names[k]: v for k, v in seqs.items()}


def test_label_counts_finds_genomes_copies_and_twins():
    labels = ["G000000001", "G000000002", "G000000002_2", "RS_GCF_000000001.1", "x"]
    counts = label_counts(labels)
    assert counts["leaves"] == 5
    assert counts["G"] == 2
    assert counts["G_like_not_matching"] == [("G#_2", 1)]
    assert counts["G_with_RS_GB_twin"] == 1
    assert genome_leaves(labels) == {"G000000001", "G000000002"}


def test_capacity_counts_only_taxa_with_enough_genomes():
    tree, _ = wol_fixture()
    labels = [leaf.label for leaf in tree.traverse_leaves()]
    lineages = lineages_by_cut(tree, (60, 30))
    genomes = genome_leaves(labels)
    capacity = clade_capacity(tree, lineages, genomes, ranks={0: 1, 1: 1})
    assert capacity[0][0] >= 1
    assert all(q <= 25 * n for n, q in capacity.values())
    assert lineage_depths(genomes, lineages) == [(2, len(genomes))]
    none = clade_capacity(tree, lineages, frozenset(labels[:3]), ranks={0: 1})
    assert none == {0: (0, 0)}


def test_dry_run_capacity_reports_disjoint_draws_and_failures():
    tree, _ = wol_fixture()
    labels = [leaf.label for leaf in tree.traverse_leaves()]
    lineages = lineages_by_cut(tree, (45, 30))
    genomes = genome_leaves(labels)
    ok = dry_run_capacity(tree, lineages, genomes, 20, ranks={0: 1})
    assert ok["clades_ok"] and ok["clade_queries"] > 0
    assert ok["singletons_wanted"] == 20 - ok["clade_queries"]
    assert ok["g_outside_clades"] < len(genomes)
    assert ok["singletons_ok"]
    too_many = dry_run_capacity(tree, lineages, genomes, 10_000, ranks={0: 1})
    assert not too_many["singletons_ok"]
    assert not dry_run_capacity(tree, lineages, genomes, 20, ranks={0: 9})["clades_ok"]


def test_length_quantiles_and_neighbours():
    _, seqs = wol_fixture()
    genomes = genome_leaves(seqs)
    rows = {g: seqs[g] for g in genomes}
    q = length_quantiles(rows)
    assert q["min"] <= q["median"] <= q["max"]
    hit, sampled = neighbour_leaves(seqs, genomes, 10)
    assert sampled == 10 and 0 <= hit <= 10
    first = sorted(genomes)[0]
    twin = {**seqs, "G999999999": seqs[first]}
    hit, sampled = neighbour_leaves(twin, frozenset({first, "G999999999"}), 1)
    assert (hit, sampled) == (1, 1)


def test_region_report_counts_hits_and_applies_g5():
    rows = {f"G{i:09d}": aligned_row("V4-150") for i in range(31)}
    rows["G000000099"] = "-" * 400
    report = region_report(rows, "V4-150")
    assert report["primer_hits"] == {"ok": 31, "no forward primer": 1}
    assert report["c_rows_kept"] == 31
    assert not report["c_dropped_g5"]
    few = region_report(dict(list(rows.items())[:5]), "V4-150")
    assert few["c_dropped_g5"]
    none = region_report({"G000000001": "-" * 400}, "V4-150")
    assert none["window"] is None and none["c_rows_kept"] == 0


def probe_inputs(capacity_ok=True):
    counts = {"G": 10_000, "G_like_not_matching": []}
    capacity = {
        "clades_ok": True, "clade_queries": 800, "singletons_ok": capacity_ok,
    }  # fmt: skip
    regions = {
        "V4-150": {"primer_hits": {"ok": 8_000}, "c_dropped_g5": False},
        "V1-V2": {"primer_hits": {"ok": 3}, "c_dropped_g5": True},
    }
    return counts, capacity, regions


def test_verdicts_pass_and_stop_on_the_preregistered_rules():
    counts, capacity, regions = probe_inputs()
    result = verdicts(counts, capacity, regions, 3000)
    assert result["G1"] == "GO N=3000" and result["G2"] == "GO"
    assert result["G3"].startswith("GO") and result["G4"].startswith("GO")
    assert "V1-V2" in result["G5"]

    counts, capacity, regions = probe_inputs(capacity_ok=False)
    counts["G"] = 5999
    counts["G_like_not_matching"] = [("G#_2", 4)]
    regions["V4-150"]["primer_hits"]["ok"] = 3000
    result = verdicts(counts, capacity, regions, 3000)
    assert result["G1"].startswith("STOP") and result["G2"].startswith("REVIEW")
    assert result["G3"].startswith("STOP") and result["G4"].startswith("STOP")


def test_probe_script_prints_every_probe(tmp_path):
    tree, seqs = wol_fixture()
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "tree.nwk").write_text(tree.newick())
    (reference / "aligned-dna-sequences.fasta").write_text(
        "".join(f">{n}\n{s}\n" for n, s in seqs.items())
    )
    lineages = lineages_by_cut(tree, (60, 30))
    taxonomy = tmp_path / "taxonomy.tsv"
    taxonomy.write_text(
        "Feature ID\tTaxon\n"
        + "".join(f"{n}\t{'; '.join(r)}\n" for n, r in lineages.items())
    )
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--reference-dir", str(reference),
            "--taxonomy", str(taxonomy), "--n-queries", "20",
            "--neighbour-sample", "5",
        ],
        capture_output=True, text=True,
        env={"PYTHONPATH": str(ROOT / "src"), "PATH": ""},
    )  # fmt: skip
    assert result.returncode == 0, result.stderr[-2000:]
    for tag in ("P1-P3", "P5", "P6", "P7", "P8", "P9", "P10", "G1", "G4", "G5"):
        assert tag in result.stdout
    assert "STOP (G1)" in result.stdout
