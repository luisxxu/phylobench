SIM10K_BACKBONE_FILES = [
    "tree.nwk",
    "aligned-dna-sequences.fasta",
    "raxml-info.txt",
    "reference_tree.nwk",
    "truth.tsv",
]


# make_sim_panel.py writes one flat directory; the queries move up a level to where the
# shared rules expect data/{dataset}/queries*.fasta.
rule sim10k_panel:
    output:
        backbone=expand("data/sim10k/B10k/{name}", name=SIM10K_BACKBONE_FILES),
        queries="data/sim10k/queries.fasta",
        aligned="data/sim10k/queries.aligned.fasta",
    params:
        python=sys.executable,
        script=SCRIPTS / "make_sim_panel.py",
        backbone_dir=lambda w, output: Path(output.backbone[0]).parent,
    log:
        "logs/sim10k/panel.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        {params.python} {params.script} {params.backbone_dir} &> {log}
        mv {params.backbone_dir}/queries.fasta {output.queries}
        mv {params.backbone_dir}/queries.aligned.fasta {output.aligned}
        """
