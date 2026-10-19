SEPP_PINNED_VERSION = "4.5.6"
SEPP_PINNED_ENV_LOCK = "envs/sepp.linux-64.pin.txt"


# run_sepp.py aligns and places in one call, so align+place is one measured phase.
# SEPP writes jplace v1; sepp_jplace.py (T26) normalizes it to v3, unmeasured.
# Scratch is cell-local (-d/-p): a tmpfs TMPDIR would count toward cgroup memory.peak.
rule sepp_pinned_place:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        info="data/{dataset}/{backbone}/raxml-info.txt",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/regions/{region}/queries.fasta",
        env_lock=SEPP_PINNED_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/raw/sepp-pinned/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/raw/sepp-pinned/align+place.run.json",
    params:
        measure=MEASURE,
        version=SEPP_PINNED_VERSION,
        python=sys.executable,
        script=SCRIPTS / "sepp_jplace.py",
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/raw/sepp-pinned/align+place.log",
    conda:
        "../../envs/sepp.yaml"
    shell:
        """
        work=$(dirname {output.jplace})/work
        rm -rf "$work" && mkdir -p "$work/out" "$work/tmp"
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            run_sepp.py -t {input.tree} -r {input.info} -a {input.ref} -f {input.queries} \
            -x {threads} -A 1000 -P 5000 -o sepp -d "$work/out" -p "$work/tmp" &> {log}
        if grep -q '"status": "ok"' {output.run_json}; then
            {params.python} {params.script} "$work/out/sepp_placement.json" {output.jplace} &>> {log}
        else
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        rm -rf "$work"
        """
