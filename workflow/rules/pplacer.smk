PPLACER_VERSION = "1.1.alpha23"
PPLACER_ENV_LOCK = "envs/pplacer.linux-64.pin.txt"


rule pplacer_place:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        info="data/{dataset}/{backbone}/raxml-info.txt",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta",
        env_lock=PPLACER_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/{aln}/pplacer/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/pplacer/place.run.json",
    params:
        measure=MEASURE,
        version=PPLACER_VERSION,
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/pplacer/place.log",
    conda:
        "../../envs/pplacer.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} --timeout 86400 \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            pplacer -j {threads} -t {input.tree} -s {input.info} -r {input.ref} \
            -o {output.jplace} {input.queries} &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        """
