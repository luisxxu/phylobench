import sys

# Private until T8 defines the shared harness constant in the Snakefile.
REFPREP_MEASURE = f"{sys.executable} -m phylobench.measure"


rule refprep_epa_ng_model:
    input:
        info="data/backbones/{backbone}/raxml-info.txt",
        aln="data/backbones/{backbone}/aligned-dna-sequences.fasta",
    output:
        "results/ref-prep/{backbone}/epa-ng.model",
    log:
        "logs/ref-prep/{backbone}/epa-ng-model.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "PYTHONPATH=src python workflow/scripts/write_epa_ng_model.py"
        " {input.info} {input.aln} {output} 2> {log}"


rule refprep_me_tree:
    input:
        tree="data/backbones/{backbone}/tree.nwk",
        aln="data/backbones/{backbone}/aligned-dna-sequences.fasta",
    output:
        tree="results/ref-prep/{backbone}/me-tree.nwk",
        run_json="results/ref-prep/{backbone}/me-tree.run.json",
    log:
        "logs/ref-prep/{backbone}/me-tree.log",
    conda:
        "../../envs/fasttree.yaml"
    threads: config["threads"]
    params:
        measure=REFPREP_MEASURE,
    shell:
        """
        OMP_NUM_THREADS={threads} {params.measure} --out {output.run_json} \
            --threads {threads} --tool-version 'FastTree 2.2.0' \
            --env-lock envs/fasttree.linux-64.pin.txt \
            -- FastTreeMP -nosupport -nt -nome -noml -intree {input.tree} {input.aln} \
            > {output.tree} 2> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            for f in {output.tree}; do [ -e "$f" ] || : > "$f"; done
        fi
        """
