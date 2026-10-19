DIPPER_VERSION = "0.1.5"
DIPPER_ENV_LOCK = "envs/dipper.linux-64.pin.txt"

localrules:
    dipper_gpu_sif,


rule dipper_m_input:
    input:
        backbone="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/{backbone}/{region}/{aln}.aligned.fasta",
    output:
        temp("results/{dataset}/{backbone}/{region}/{aln}/dipper-m.in.fasta"),
    wildcard_constraints:
        aln="hmmalign|oracle",
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/dipper-m.in.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        # awk 1 adds the final newline a file may lack, so records never join.
        "awk 1 {input.backbone} {input.queries} > {output} 2> {log}"


rule dipper_r_input:
    input:
        backbone="data/{dataset}/{backbone}/aligned-dna-sequences.fasta",
        queries="results/{dataset}/regions/{region}/queries.fasta",
        script=SCRIPTS / "ungap_concat.awk",
    output:
        temp("results/{dataset}/{backbone}/{region}/raw/dipper-r.in.fasta"),
    log:
        "results/{dataset}/{backbone}/{region}/raw/dipper-r.in.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "awk -f {input.script} {input.backbone} {input.queries} > {output} 2> {log}"


rule dipper_m_place:
    input:
        fasta="results/{dataset}/{backbone}/{region}/{aln}/dipper-m.in.fasta",
        tree="data/{dataset}/{backbone}/tree.nwk",
        env_lock=DIPPER_ENV_LOCK,
    output:
        nwk="results/{dataset}/{backbone}/{region}/{aln}/dipper-m/dipper.nwk",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/dipper-m/place.run.json",
    wildcard_constraints:
        aln="hmmalign|oracle",
    params:
        measure=MEASURE,
        version=DIPPER_VERSION,
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/dipper-m/place.log",
    conda:
        "../../envs/dipper.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            dipper_cpu --add -i m -d 2 -T {threads} -t {input.tree} \
            -I {input.fasta} -O {output.nwk} &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.nwk} ] || : > {output.nwk}
        fi
        """


rule dipper_r_place:
    input:
        fasta="results/{dataset}/{backbone}/{region}/raw/dipper-r.in.fasta",
        tree="data/{dataset}/{backbone}/tree.nwk",
        env_lock=DIPPER_ENV_LOCK,
    output:
        nwk="results/{dataset}/{backbone}/{region}/raw/dipper-r/dipper.nwk",
        run_json="results/{dataset}/{backbone}/{region}/raw/dipper-r/place.run.json",
    params:
        measure=MEASURE,
        version=DIPPER_VERSION,
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/raw/dipper-r/place.log",
    conda:
        "../../envs/dipper.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --threads {threads} \
            --tool-version '{params.version}' --env-lock {input.env_lock} -- \
            dipper_cpu --add -i r -T {threads} -t {input.tree} \
            -I {input.fasta} -O {output.nwk} &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.nwk} ] || : > {output.nwk}
        fi
        """


rule dipper_gpu_sif:
    input:
        "envs/dipper-gpu.image.txt",
    output:
        "data/images/dipper-gpu.sif",
    log:
        "data/images/dipper-gpu.pull.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        'apptainer pull {output} "$(cat {input})" &> {log}'


rule dipper_m_gpu_place:
    input:
        fasta="results/{dataset}/{backbone}/{region}/{aln}/dipper-m.in.fasta",
        tree="data/{dataset}/{backbone}/tree.nwk",
        sif="data/images/dipper-gpu.sif",
        image="envs/dipper-gpu.image.txt",
    output:
        nwk="results/{dataset}/{backbone}/{region}/{aln}/dipper-m-gpu/dipper.nwk",
        run_json="results/{dataset}/{backbone}/{region}/{aln}/dipper-m-gpu/place.run.json",
    wildcard_constraints:
        aln="hmmalign|oracle",
    params:
        measure=MEASURE,
    threads: config["threads"]
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/dipper-m-gpu/place.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        {params.measure} --gpu --out {output.run_json} --threads {threads} \
            --tool-version "$(cat {input.image})" --env-lock {input.image} -- \
            apptainer exec --nv {input.sif} dipper --add -i m -d 2 \
            -t {input.tree} -I {input.fasta} -O {output.nwk} &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.nwk} ] || : > {output.nwk}
        fi
        """


rule dipper_to_jplace:
    input:
        gate="results/{dataset}/{backbone}/{region}/{aln}/{arm}/place.run.json",
        nwk="results/{dataset}/{backbone}/{region}/{aln}/{arm}/dipper.nwk",
        backbone="data/{dataset}/{backbone}/tree.nwk",
        region_queries="results/{dataset}/regions/{region}/queries.fasta",
    output:
        "results/{dataset}/{backbone}/{region}/{aln}/{arm}/placements.jplace",
    wildcard_constraints:
        arm="dipper-m|dipper-m-gpu|dipper-r",
    params:
        python=sys.executable,
        script=SCRIPTS / "dipper_to_jplace.py",
    threads: 1
    log:
        "results/{dataset}/{backbone}/{region}/{aln}/{arm}/to_jplace.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --gate-run-json {input.gate} "
        "--dipper-tree {input.nwk} --backbone-tree {input.backbone} "
        "--region-queries {input.region_queries} --out {output} &> {log}"
