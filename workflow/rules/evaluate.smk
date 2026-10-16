from importlib.metadata import version as package_version
import json
import shlex

from phylobench.tables import run_specs

PHYLOBENCH_VERSION = f"phylobench {package_version('phylobench')}"


wildcard_constraints:
    aln="[^/]+",
    arm="[^/]+",


rule similarity:
    input:
        queries="results/{dataset}/{backbone}/{region}/hmmalign.aligned.fasta",
        backbone="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
    output:
        "results/{dataset}/{backbone}/{region}/similarity.tsv",
    params:
        python=sys.executable,
        script=SCRIPTS / "similarity.py",
    log:
        "results/{dataset}/{backbone}/{region}/similarity.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --queries {input.queries} "
        "--backbone {input.backbone} --out {output} &> {log}"


rule evaluate:
    input:
        gate=lambda wildcards: (
            f"results/{wildcards.dataset}/{wildcards.backbone}/{wildcards.region}/"
            f"{wildcards.aln}/{wildcards.arm}/"
            f"{next(iter(config['arms'][wildcards.arm]['phases'].values()))}"
        ),
        jplace="results/{dataset}/{backbone}/{region}/{aln}/{arm}/placements.jplace",
        backbone="data/{dataset}/{backbone}/tree.nwk",
        reference="data/{dataset}/{backbone}/reference_tree.nwk",
        region_queries="results/{dataset}/regions/{region}/queries.fasta",
        similarity="results/{dataset}/{backbone}/{region}/similarity.tsv",
    output:
        accuracy="results/{dataset}/{backbone}/{region}/{aln}/{arm}/accuracy.tsv",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/{arm}/evaluate.run.json",
    params:
        measure=MEASURE,
        version=PHYLOBENCH_VERSION,
        python=sys.executable,
        script=SCRIPTS / "evaluate_run.py",
        run_id="{dataset}/{backbone}/{region}/{aln}/{arm}",
    threads: 1
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/{arm}/evaluate.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock envs/phylobench.yaml -- \
            {params.python} {params.script} --run-id {params.run_id} \
            --gate-run-json {input.gate} --jplace {input.jplace} \
            --backbone-tree {input.backbone} --reference-tree {input.reference} \
            --region-queries {input.region_queries} --similarity {input.similarity} \
            --out {output.accuracy} &> {log}
        grep -q '"status": "ok"' {output.run_json}
        """


rule collect_tables:
    input:
        lambda wildcards: [
            path
            for s in run_specs(config, config["arms"])
            for path in (
                [s.run_json, str(Path(s.run_json).with_name("accuracy.tsv"))]
                if s.phase == "evaluate"
                else [s.run_json]
            )
        ],
    output:
        runs="results/runs.tsv",
        accuracy="results/accuracy.tsv",
    params:
        python=sys.executable,
        script=SCRIPTS / "collect_tables.py",
        specs=lambda wildcards: shlex.quote(
            json.dumps([s._asdict() for s in run_specs(config, config["arms"])])
        ),
    log:
        "results/collect_tables.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --specs {params.specs} "
        "--out-runs {output.runs} --out-accuracy {output.accuracy} &> {log}"
