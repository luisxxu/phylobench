UPP_VERSION = "sepp 4.5.6 (run_upp.py)"
UPP_ENV_LOCK = "envs/sepp.linux-64.pin.txt"
# 1440 min (the prep-rule runtime in profiles/slurm-pilot) minus 40 min for the harness to
# write its DNF row; tests/test_profile.py ties this to both profiles.
UPP_TIMEOUT_S = 84000


# UPP aligns the queries to the fixed backbone alignment and tree (-a, -t), so PASTA never
# runs. -d/-p must be absolute. Scratch is cell-local: a tmpfs TMPDIR would count toward
# cgroup memory.peak.
rule upp_align:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/regions/{region}/queries.fasta",
        env_lock=UPP_ENV_LOCK,
    output:
        masked="results/{dataset}/{backbone}/{region}/align/upp.masked.fasta",
        insertion_columns="results/{dataset}/{backbone}/{region}/align/upp.insertion_columns.txt",
        run_json="results/{dataset}/{backbone}/{region}/align/upp.run.json",
    params:
        measure=MEASURE,
        version=UPP_VERSION,
        timeout=UPP_TIMEOUT_S,
        workdir=lambda wildcards, output: str(Path(output.run_json).parent / "upp"),
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/align/upp.log",
    conda:
        "../../envs/sepp.yaml"
    shell:
        """
        work="$PWD/{params.workdir}"
        rm -rf "$work" && mkdir -p "$work/out" "$work/tmp"
        {params.measure} --out {output.run_json} --threads {threads} --timeout {params.timeout} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            run_upp.py -s {input.queries} -a {input.ref} -t {input.tree} -A 100 -m dna \
            -x {threads} -o upp -d "$work/out" -p "$work/tmp" &> {log}
        if grep -q '"status": "ok"' {output.run_json}; then
            test -s "$work/out/upp_alignment_masked.fasta"
            mv "$work/out/upp_alignment_masked.fasta" {output.masked}
            mv "$work/out/upp_insertion_columns.txt" {output.insertion_columns}
        else
            # keeps the DAG going so the DNF row is recorded, as mask_hmmalign does
            : > {output.masked}
            : > {output.insertion_columns}
        fi
        rm -rf "$work/tmp"
        """


rule upp_queries:
    input:
        masked="results/{dataset}/{backbone}/{region}/align/upp.masked.fasta",
        run_json="results/{dataset}/{backbone}/{region}/align/upp.run.json",
        backbone="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/regions/{region}/queries.fasta",
    output:
        "results/{dataset}/{backbone}/{region}/upp.aligned.fasta",
    params:
        python=sys.executable,
        script=SCRIPTS / "upp_queries.py",
    log:
        "results/{dataset}/{backbone}/{region}/align/upp_queries.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --masked {input.masked} "
        "--backbone {input.backbone} --queries {input.queries} "
        "--run-json {input.run_json} --out {output} &> {log}"
