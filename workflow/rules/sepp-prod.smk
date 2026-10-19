SEPP_PROD_VERSION = "sepp 4.5.6, pplacer 1.1.alpha19"
SEPP_PROD_ENV_LOCK = "envs/sepp-prod.linux-64.pin.txt"


# run-sepp.sh overrides -d and -p, so the per-job work dir is the cwd and TMPDIR is
# set explicitly. The jplace converter runs outside the harness (unmeasured).
rule sepp_prod:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        info="data/{dataset}/{backbone}/raxml-info.txt",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/regions/{region}/queries.fasta",
        env_lock=SEPP_PROD_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/raw/sepp-prod/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/raw/sepp-prod/sepp.run.json",
    params:
        measure=MEASURE,
        version=SEPP_PROD_VERSION,
        python=sys.executable,
        script=SCRIPTS / "sepp_jplace.py",
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/raw/sepp-prod/sepp.log",
    conda:
        "../../envs/sepp-prod.yaml"
    shell:
        """
        root=$PWD; work=$root/$(dirname {output.jplace})/work
        rm -rf $work; mkdir -p $work/tmp; cd $work
        export TMPDIR=$work/tmp
        {params.measure} --out $root/{output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock $root/{input.env_lock} -- \
            run-sepp.sh $root/{input.queries} sepp -x {threads} -A 1000 -P 5000 \
            -a $root/{input.ref} -t $root/{input.tree} -r $root/{input.info} &> $root/{log}
        cat sepp-sepp-*.log >> $root/{log} 2>/dev/null || true
        cd $root
        if grep -q '"status": "ok"' {output.run_json}; then
            {params.python} {params.script} $work/sepp_placement.json {output.jplace} &>> {log}
        else
            : > {output.jplace}
        fi
        """
