HMMER_VERSION = "hmmer 3.4"
HMMER_ENV_LOCK = "envs/hmmer.linux-64.pin.txt"
SCRIPTS = Path(workflow.basedir) / "scripts"


wildcard_constraints:
    dataset="[^/]+",
    backbone="[^/]+",
    region="full-length|V4-150|V3-V4|V1-V2|V1-V3|V4-V5|V5-V6|V6-V8|V7-V9",


# Python rules call the phylobench env's interpreter by absolute path (as MEASURE does),
# so the package resolves under --use-conda too.
rule extract_regions:
    input:
        queries="data/{dataset}/queries.fasta",
    output:
        fasta="results/{dataset}/regions/{region}/queries.fasta",
        exclusions="results/{dataset}/regions/{region}/exclusions.tsv",
        counts="results/{dataset}/regions/{region}/exclusion_counts.tsv",
    params:
        python=sys.executable,
        script=SCRIPTS / "extract_regions.py",
    log:
        "results/{dataset}/regions/{region}/extract_regions.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --region {wildcards.region} "
        "--queries {input.queries} --out-fasta {output.fasta} "
        "--out-exclusions {output.exclusions} --out-counts {output.counts} &> {log}"


rule oracle_region:
    input:
        queries="data/{dataset}/queries.fasta",
        aligned="data/{dataset}/queries.aligned.fasta",
        backbone="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
    output:
        "results/{dataset}/{backbone}/{region}/oracle.aligned.fasta",
    params:
        python=sys.executable,
        script=SCRIPTS / "oracle_region.py",
    log:
        "results/{dataset}/{backbone}/{region}/oracle_region.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --region {wildcards.region} "
        "--queries {input.queries} --aligned {input.aligned} "
        "--backbone {input.backbone} --out {output} &> {log}"


rule backbone_stockholm:
    input:
        "data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
    output:
        "results/{dataset}/{backbone}/align/backbone.sto",
    params:
        python=sys.executable,
        script=SCRIPTS / "backbone_stockholm.py",
    log:
        "results/{dataset}/{backbone}/align/backbone_stockholm.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --aligned {input} --out {output} &> {log}"


rule hmmbuild:
    input:
        sto="results/{dataset}/{backbone}/align/backbone.sto",
        env_lock=HMMER_ENV_LOCK,
    output:
        hmm="results/{dataset}/{backbone}/align/backbone.hmm",
        run_json="results/{dataset}/{backbone}/align/hmmbuild.run.json",
    params:
        measure=MEASURE,
        version=HMMER_VERSION,
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/align/hmmbuild.log",
    conda:
        "../../envs/hmmer.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            hmmbuild --dna --hand --cpu {threads} {output.hmm} {input.sto} &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.hmm} ] || : > {output.hmm}
        fi
        """


# hmmalign 3.4 has no thread flag (unlike hmmbuild); the 16 threads are requested
# from the scheduler and recorded in run.json so the cost comparison stays honest.
rule hmmalign:
    input:
        hmm="results/{dataset}/{backbone}/align/backbone.hmm",
        queries="results/{dataset}/regions/{region}/queries.fasta",
        env_lock=HMMER_ENV_LOCK,
    output:
        afa="results/{dataset}/{backbone}/{region}/align/hmmalign.afa",
        run_json="results/{dataset}/{backbone}/{region}/align/hmmalign.run.json",
    params:
        measure=MEASURE,
        version=HMMER_VERSION,
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/align/hmmalign.log",
    conda:
        "../../envs/hmmer.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            hmmalign --dna --outformat afa -o {output.afa} {input.hmm} {input.queries} &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.afa} ] || : > {output.afa}
        fi
        """


rule mask_hmmalign:
    input:
        afa="results/{dataset}/{backbone}/{region}/align/hmmalign.afa",
        run_json="results/{dataset}/{backbone}/{region}/align/hmmalign.run.json",
        backbone="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
    output:
        "results/{dataset}/{backbone}/{region}/hmmalign.aligned.fasta",
    params:
        python=sys.executable,
        script=SCRIPTS / "mask_hmmalign.py",
    log:
        "results/{dataset}/{backbone}/{region}/align/mask_hmmalign.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --afa {input.afa} "
        "--backbone {input.backbone} --run-json {input.run_json} --out {output} &> {log}"
