import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
BACKBONE_FILES = ("tree.nwk", "aligned-dna-sequences.fasta", "raxml-info.txt")
DATASETS = ("emp-10k", "emp-100k")
BACKBONES = ("B10k", "Bfull")
REGION = "V4-150"


def run_script(script, *args, cwd=None):
    return subprocess.run(
        [sys.executable, str(SCRIPTS / script), *map(str, args)],
        capture_output=True,
        text=True,
        cwd=cwd,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def write_gg2(root: Path, n_per_bin: int = 2) -> None:
    """data/gg2 backbone files and the MANIFEST.tsv build_gg2.py writes for them."""
    rows = []
    for backbone in BACKBONES:
        for name in BACKBONE_FILES:
            path = write(root / "data" / "gg2" / backbone / name, f"{backbone}{name}\n")
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            n = "" if name == "raxml-info.txt" else 7
            rows.append(
                ("file", backbone, "", "", f"data/gg2/{backbone}/{name}", n, sha, "")
            )
    columns = "kind backbone region bin path n sha256 flagged".split()
    pd.DataFrame(rows, columns=columns).to_csv(
        root / "data" / "gg2" / "MANIFEST.tsv", sep="\t", index=False
    )


def link_emp_backbones(root: Path) -> None:
    for dataset in DATASETS:
        for backbone in BACKBONES:
            for name in BACKBONE_FILES:
                link = root / "data" / dataset / backbone / name
                link.parent.mkdir(parents=True, exist_ok=True)
                link.symlink_to(Path("../../gg2") / backbone / name)


def write_biom(path: Path, ids: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.attrs["format-version"] = np.array([2, 1])
        f.create_dataset("observation/ids", data=np.array(ids, dtype="S"))


def write_emp_stage(root: Path) -> None:
    """Everything emp_manifest.py reads, for both datasets and both backbones."""
    write_gg2(root)
    link_emp_backbones(root)
    write_biom(
        root / "data" / "emp" / "emp_deblur_150bp.release1.biom", ["A" * 150] * 3
    )
    for dataset in DATASETS:
        write(
            root / "results" / dataset / "regions" / REGION / "queries.fasta",
            ">asv000001\nACGT\n>asv000002\nACGT\n",
        )
        for backbone in BACKBONES:
            sims = pd.DataFrame(
                {"query": ["asv000001", "asv000002"], "sim": 0.9, "bin": "90-95"}
            )
            sims_dir = root / "results" / dataset / backbone / REGION
            sims_dir.mkdir(parents=True, exist_ok=True)
            sims.to_csv(sims_dir / "similarity.tsv", sep="\t", index=False)


def run_manifest(root: Path):
    return run_script(
        "emp_manifest.py",
        "--datasets", *DATASETS,
        "--backbones", *BACKBONES,
        "--out", root / "data" / "emp" / "MANIFEST.tsv",
        cwd=root,
    )  # fmt: skip


def write_run_json(path: Path, status: str = "ok") -> Path:
    return write(path, json.dumps({"status": status}))
