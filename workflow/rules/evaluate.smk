from importlib.metadata import version as package_version
import json
from pathlib import Path

from phylobench.tables import run_specs, score_table

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


rule concordance:
    input:
        gate=lambda wildcards: (
            f"results/{wildcards.dataset}/{wildcards.backbone}/{wildcards.region}/"
            f"{wildcards.aln}/{wildcards.arm}/"
            f"{next(iter(config['arms'][wildcards.arm]['phases'].values()))}"
        ),
        jplace="results/{dataset}/{backbone}/{region}/{aln}/{arm}/placements.jplace",
        # The sepp_prod rule's outputs, spelled out: config['arms'] lacks sepp-prod under
        # --config arm=X, so the arm gate lambda cannot name them.
        reference_gate="results/{dataset}/{backbone}/{region}/raw/sepp-prod/sepp.run.json",
        reference_jplace="results/{dataset}/{backbone}/{region}/raw/sepp-prod/placements.jplace",
        backbone="data/{dataset}/{backbone}/tree.nwk",
        region_queries="results/{dataset}/regions/{region}/queries.fasta",
        similarity="results/{dataset}/{backbone}/{region}/similarity.tsv",
    output:
        concordance="results/{dataset}/{backbone}/{region}/{aln}/{arm}/concordance.tsv",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/{arm}/concordance.run.json",
    wildcard_constraints:
        dataset="emp-10k|emp-100k",
    params:
        measure=MEASURE,
        version=PHYLOBENCH_VERSION,
        python=sys.executable,
        script=SCRIPTS / "concordance_run.py",
        run_id="{dataset}/{backbone}/{region}/{aln}/{arm}",
        reference_run_id="{dataset}/{backbone}/{region}/raw/sepp-prod",
    threads: 1
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/{arm}/concordance.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock envs/phylobench.yaml -- \
            {params.python} {params.script} --run-id {params.run_id} \
            --reference-run-id {params.reference_run_id} \
            --gate-run-json {input.gate} --reference-gate-run-json {input.reference_gate} \
            --jplace {input.jplace} --reference-jplace {input.reference_jplace} \
            --backbone-tree {input.backbone} --region-queries {input.region_queries} \
            --similarity {input.similarity} --out {output.concordance} &> {log}
        grep -q '"status": "ok"' {output.run_json}
        """


# The specs travel in params so an exclude or config change reruns this rule; a single
# argv string holding them would exceed Linux's 128 KiB MAX_ARG_STRLEN at full scale.
localrules:
    run_specs,


rule run_specs:
    output:
        "results/run_specs.json",
    params:
        specs=lambda wildcards: [s._asdict() for s in run_specs(config, config["arms"])],
    log:
        "results/run_specs.log",
    run:
        Path(output[0]).write_text(json.dumps(params.specs))


rule collect_tables:
    input:
        specs="results/run_specs.json",
        run_jsons=lambda wildcards: [
            path
            for s in run_specs(config, config["arms"])
            for path in (
                [s.run_json, score_table(s)] if s.phase == "evaluate" else [s.run_json]
            )
        ],
    output:
        runs="results/runs.tsv",
        accuracy="results/accuracy.tsv",
        concordance="results/concordance.tsv",
    params:
        python=sys.executable,
        script=SCRIPTS / "collect_tables.py",
    log:
        "results/collect_tables.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --specs {input.specs} "
        "--out-runs {output.runs} --out-accuracy {output.accuracy} "
        "--out-concordance {output.concordance} &> {log}"
