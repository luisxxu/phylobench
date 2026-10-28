localrules:
    plots,


FIGURES = [
    "accuracy_edge_dist",
    "accuracy_exact",
    "accuracy_delta_err",
    "verdict",
    "oracle_gap",
    "cost_vs_similarity",
    "cost_vs_backbone",
    "cost_vs_query_count",
    "end_to_end",
    "gpu_memory",
]
ACCURACY_FIGURES = FIGURES[:5]
COST_FIGURES = FIGURES[5:]
TABLES = config["plot_tables"]


# Not in `all`: real tables are plotted by running the scripts directly (README).
rule plots:
    input:
        [f"results/figures/{name}.png" for name in FIGURES],


rule plot_accuracy:
    input:
        accuracy=f"{TABLES}/accuracy.tsv",
    output:
        [f"results/figures/{name}.png" for name in ACCURACY_FIGURES],
    params:
        python=sys.executable,
        script=SCRIPTS / "plot_accuracy.py",
        outdir=lambda wildcards, output: Path(output[0]).parent,
    log:
        "results/figures/plot_accuracy.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --accuracy {input.accuracy} "
        "--outdir {params.outdir} &> {log}"


rule plot_cost:
    input:
        runs=f"{TABLES}/runs.tsv",
        accuracy=f"{TABLES}/accuracy.tsv",
    output:
        [f"results/figures/{name}.png" for name in COST_FIGURES],
    params:
        python=sys.executable,
        script=SCRIPTS / "plot_cost.py",
        outdir=lambda wildcards, output: Path(output[0]).parent,
    log:
        "results/figures/plot_cost.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --runs {input.runs} "
        "--accuracy {input.accuracy} --outdir {params.outdir} &> {log}"
