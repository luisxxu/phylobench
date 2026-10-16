# The fixture's backbone name is fixed by config/smoke.yaml (backbones: [B100]).
rule stage_tiny_backbone:
    input:
        "tests/fixtures/tiny/{name}",
    output:
        "data/tiny/B100/{name}",
    wildcard_constraints:
        name="tree.nwk|aligned-dna-sequences.fasta|raxml-info.txt|reference_tree.nwk",
    log:
        "logs/stage_tiny/{name}.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "cp {input} {output} &> {log}"


rule stage_tiny_queries:
    input:
        "tests/fixtures/tiny/{name}",
    output:
        "data/tiny/{name}",
    wildcard_constraints:
        name="queries.fasta|queries.aligned.fasta",
    log:
        "logs/stage_tiny/{name}.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "cp {input} {output} &> {log}"


# Proves a non-ok run becomes a DNF row (docs/harness.md section 7).
rule forced_dnf:
    output:
        jplace="results/tiny/B100/full-length/hmmalign/forced-dnf/placements.jplace",
        run_json="results/tiny/B100/full-length/hmmalign/forced-dnf/place.run.json",
    params:
        measure=MEASURE,
    log:
        "results/tiny/B100/full-length/hmmalign/forced-dnf/place.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        {params.measure} --out {output.run_json} --timeout 1 --threads 1 -- \
            sleep 30 &> {log}
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.jplace} ] || : > {output.jplace}
        fi
        """
