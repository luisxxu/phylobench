EPA_NG_ENV_LOCK = "envs/epa-ng.linux-64.pin.txt"


rule epa_ng_place:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        model="results/ref-prep/{dataset}/{backbone}/epa-ng.model",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta",
        env_lock=EPA_NG_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/{aln}/epa-ng/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/epa-ng/place.run.json",
    params:
        measure=MEASURE,
        version="epa-ng 0.3.8",
        outdir=lambda wildcards, output: str(Path(output.jplace).parent),
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/epa-ng/place.log",
    conda:
        "../../envs/epa-ng.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            epa-ng --ref-msa {input.ref} --tree {input.tree} --query {input.queries} \
            --model {input.model} --outdir {params.outdir} --threads {threads} --redo &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        if grep -q '"status": "ok"' {output.run_json}; then
            mv {params.outdir}/epa_result.jplace {output.jplace}
        fi
        """
