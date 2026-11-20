import argparse
from functools import cache
from pathlib import Path

import pandas as pd
import yaml

from phylobench.accuracy import prepare_backbone, true_positions
from phylobench.jplace import parse_jplace
from phylobench.regions import read_fasta
from phylobench.skeleton import skeleton_map, skeleton_newick, skeleton_scores
from phylobench.tables import load_arms, plan_cells
from phylobench.wol import genome_leaves

ROOT = Path(__file__).parents[2]


@cache
def load_backbone(dataset: str, backbone: str):
    newick = Path(f"data/{dataset}/{backbone}/tree.nwk").read_text()
    full = prepare_backbone(newick)
    skeleton = prepare_backbone(skeleton_newick(newick, genome_leaves(full.hashes)))
    return full, skeleton, skeleton_map(newick, skeleton)


@cache
def load_truth(dataset: str, backbone: str, region: str):
    queries = [
        name
        for name, _ in read_fasta(f"results/{dataset}/regions/{region}/queries.fasta")
    ]
    reference = Path(f"data/{dataset}/{backbone}/reference_tree.nwk").read_text()
    skeleton = load_backbone(dataset, backbone)[1]
    return queries, *true_positions(reference, queries, skeleton)


def main():
    parser = argparse.ArgumentParser(
        description="Score every cell with accuracy rows on the G-only skeleton."
    )
    parser.add_argument("configs", type=Path, nargs="+")
    parser.add_argument("--dataset", default="gg2-wol")
    parser.add_argument("--out", type=Path, default="results/skeleton.tsv")
    args = parser.parse_args()

    config = {}
    for path in args.configs:
        config.update(yaml.safe_load(path.read_text()))
    planned = [
        c
        for c in plan_cells(config, load_arms(ROOT / "config" / "arms"))
        if c.dataset == args.dataset
    ]

    rows, skipped = [], 0
    for cell in planned:
        cell_dir = Path("results", cell.dataset, cell.backbone, cell.region)
        cell_dir = cell_dir / cell.aln / cell.arm
        if len(pd.read_csv(cell_dir / "accuracy.tsv", sep="\t")) == 0:
            skipped += 1
            continue
        backbone, skeleton, node_map = load_backbone(cell.dataset, cell.backbone)
        queries, true_node, true_weight = load_truth(
            cell.dataset, cell.backbone, cell.region
        )
        doc = parse_jplace((cell_dir / "placements.jplace").read_text())
        scores = skeleton_scores(
            doc, backbone, skeleton, node_map, true_node, true_weight, queries
        )
        rows.append(pd.DataFrame(scores).assign(run_id=cell.run_id))
    table = pd.concat(rows, ignore_index=True)[["run_id", "query", "edge_dist"]]
    table.to_csv(args.out, sep="\t", index=False)
    print(f"{len(rows)} scored, {skipped} skipped (no accuracy rows)")


if __name__ == "__main__":
    main()
