import argparse
import json
from pathlib import Path

from phylobench.dipper_adapter import dipper_to_jplace
from phylobench.regions import read_fasta


def main():
    parser = argparse.ArgumentParser(
        description="Convert DIPPER's Newick output to a jplace for evaluation."
    )
    parser.add_argument("--gate-run-json", type=Path, required=True)
    parser.add_argument("--dipper-tree", type=Path, required=True)
    parser.add_argument("--backbone-tree", type=Path, required=True)
    parser.add_argument("--region-queries", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    # A DNF'd placement left an empty tree; its run.json carries the failure.
    if json.loads(args.gate_run_json.read_text())["status"] != "ok":
        args.out.write_text("")
        return
    queries = [name for name, _ in read_fasta(args.region_queries)]
    doc = dipper_to_jplace(
        args.dipper_tree.read_text(), args.backbone_tree.read_text(), queries
    )
    # NaN lengths are written as bare NaN: evaluation-only, never for guppy.
    args.out.write_text(json.dumps(doc))


if __name__ == "__main__":
    main()
