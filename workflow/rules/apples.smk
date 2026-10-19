APPLES_VERSION = "2.0.11"
APPLES_ENV_LOCK = "envs/apples.linux-64.pin.txt"


rule apples2_place:
    input:
        tree="results/ref-prep/{dataset}/{backbone}/me-tree.nwk",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta",
        env_lock=APPLES_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/{aln}/apples2/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/apples2/place.run.json",
    params:
        measure=MEASURE,
        version=APPLES_VERSION,
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/apples2/place.log",
    conda:
        "../../envs/apples.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            run_apples.py -s {input.ref} -q {input.queries} -t {input.tree} \
            -T {threads} -D -o {output.jplace} &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        """
