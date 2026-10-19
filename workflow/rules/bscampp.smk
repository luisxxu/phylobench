BSCAMPP_ENV_LOCK = "envs/bscampp.linux-64.pin.txt"


rule bscampp_place:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        model="results/ref-prep/{dataset}/{backbone}/epa-ng.model",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta",
        env_lock=BSCAMPP_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/{aln}/bscampp/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/bscampp/place.run.json",
    params:
        measure=MEASURE,
        version="bscampp 1.0.8 (epa-ng 0.3.8)",
        outdir=lambda wildcards, output: str(Path(output.jplace).parent),
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/bscampp/place.log",
    conda:
        "../../envs/bscampp.yaml"
    shell:
        """
        rm -rf {params.outdir}/tmp0
        OMP_NUM_THREADS={threads} {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            run_bscampp.py --placement-method epa-ng -i {input.model} \
            -t {input.tree} -a {input.ref} -q {input.queries} \
            -d {params.outdir} -o placements.jplace --threads {threads} \
            --cpus-per-job 2 -b 2000 -V 5 &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        if grep -q '"status": "ok"' {output.run_json}; then
            grep -qF "main.config: $CONDA_PREFIX/share/bscampp/main.config" {log}
            grep -qF "Configs.epang_path: $CONDA_PREFIX/bin/epa-ng" {log}
        fi
        if grep -Eq 'Cannot find home.path|Cannot find main.config|Initializing the config file' {log}; then
            exit 1
        fi
        """


rule epa_ng_scampp_place:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        model="results/ref-prep/{dataset}/{backbone}/epa-ng.model",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta",
        env_lock=BSCAMPP_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/{aln}/epa-ng-scampp/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/epa-ng-scampp/place.run.json",
    params:
        measure=MEASURE,
        version="bscampp 1.0.8 (epa-ng 0.3.8)",
        outdir=lambda wildcards, output: str(Path(output.jplace).parent),
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/epa-ng-scampp/place.log",
    conda:
        "../../envs/bscampp.yaml"
    shell:
        """
        rm -rf {params.outdir}/tmp0
        OMP_NUM_THREADS={threads} {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            run_scampp.py --placement-method epa-ng -i {input.model} \
            -t {input.tree} -a {input.ref} -q {input.queries} \
            -d {params.outdir} -o placements.jplace --threads {threads} \
            --cpus-per-job 1 -b 2000 &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        if grep -q '"status": "ok"' {output.run_json}; then
            grep -qF "main.config: $CONDA_PREFIX/share/bscampp/main.config" {log}
            grep -qF "Configs.epang_path: $CONDA_PREFIX/bin/epa-ng" {log}
        fi
        if grep -Eq 'Cannot find home.path|Cannot find main.config|Initializing the config file' {log}; then
            exit 1
        fi
        """


rule pplacer_scampp_place:
    input:
        tree="data/{dataset}/{backbone}/tree.nwk",
        model="data/{dataset}/{backbone}/raxml-info.txt",
        ref="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta",
        env_lock=BSCAMPP_ENV_LOCK,
    output:
        jplace="results/{dataset}/{backbone}/{region}/{aln}/pplacer-scampp/placements.jplace",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/pplacer-scampp/place.run.json",
    params:
        measure=MEASURE,
        version="bscampp 1.0.8 (pplacer 1.1.alpha23)",
        outdir=lambda wildcards, output: str(Path(output.jplace).parent),
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/pplacer-scampp/place.log",
    conda:
        "../../envs/bscampp.yaml"
    shell:
        """
        rm -rf {params.outdir}/tmp0
        OMP_NUM_THREADS={threads} {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            run_scampp.py --placement-method pplacer -i {input.model} \
            -t {input.tree} -a {input.ref} -q {input.queries} \
            -d {params.outdir} -o placements.jplace --threads {threads} \
            --cpus-per-job 1 -b 2000 &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        if grep -q '"status": "ok"' {output.run_json}; then
            grep -qF "main.config: $CONDA_PREFIX/share/bscampp/main.config" {log}
            grep -qF "Configs.pplacer_path: $CONDA_PREFIX/bin/pplacer" {log}
            grep -qF "Configs.taxit_path: $CONDA_PREFIX/bin/taxit" {log}
        fi
        if grep -Eq 'Cannot find home.path|Cannot find main.config|Initializing the config file' {log}; then
            exit 1
        fi
        """
