GG2_URL = "https://ftp.microbio.me/greengenes_release/2022.10"
GG2_REFERENCE_SHA256 = "fa5d32f936a39735a4f99552044ed4fcbae0d91fb48a679f9d329e2a708eb973"
GG2_TAXONOMY_SHA256 = "6495d6da4097ff51391ae7ce8df09c0ee00807c8bb64e4ce8cae1c9d5e07f5f0"
GG2_REFERENCE_MEMBERS = ["tree.nwk", "aligned-dna-sequences.fasta", "raxml-info.txt"]
GG2_BACKBONE_FILES = [
    "tree.nwk",
    "aligned-dna-sequences.fasta",
    "raxml-info.txt",
    "reference_tree.nwk",
]


rule gg2_download:
    output:
        reference="data/gg2/download/2022.10.backbone.sepp-reference.qza",
        taxonomy="data/gg2/download/2022.10.backbone.tax.qza",
    params:
        fetch=SCRIPTS / "fetch.sh",
        reference_url=f"{GG2_URL}/2022.10.backbone.sepp-reference.qza",
        taxonomy_url=f"{GG2_URL}/2022.10.backbone.tax.qza",
        reference_sha256=GG2_REFERENCE_SHA256,
        taxonomy_sha256=GG2_TAXONOMY_SHA256,
    log:
        "logs/gg2/download.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        {params.fetch} {params.reference_url} {params.reference_sha256} \
            {output.reference} &> {log}
        {params.fetch} {params.taxonomy_url} {params.taxonomy_sha256} \
            {output.taxonomy} &>> {log}
        """


rule gg2_unzip:
    input:
        reference="data/gg2/download/2022.10.backbone.sepp-reference.qza",
        taxonomy="data/gg2/download/2022.10.backbone.tax.qza",
    output:
        reference=expand(
            "data/gg2/shipped/sepp-reference/{member}", member=GG2_REFERENCE_MEMBERS
        ),
        taxonomy="data/gg2/shipped/tax/taxonomy.tsv",
    params:
        python=sys.executable,
        script=SCRIPTS / "unzip_qza.py",
        reference_dir=lambda w, output: Path(output.reference[0]).parent,
        taxonomy_dir=lambda w, output: Path(output.taxonomy).parent,
        members=" ".join(GG2_REFERENCE_MEMBERS),
    log:
        "logs/gg2/unzip.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        {params.python} {params.script} {input.reference} {params.reference_dir} \
            {params.members} &> {log}
        {params.python} {params.script} {input.taxonomy} {params.taxonomy_dir} \
            taxonomy.tsv &>> {log}
        """


rule gg2_build:
    input:
        reference=rules.gg2_unzip.output.reference,
        taxonomy=rules.gg2_unzip.output.taxonomy,
        downloads=rules.gg2_download.output,
    output:
        manifest="data/gg2/MANIFEST.tsv",
        backbones=expand(
            "data/gg2/{backbone}/{name}",
            backbone=["Bfull", "B10k"],
            name=GG2_BACKBONE_FILES,
        ),
        queries="data/gg2/queries.fasta",
        aligned="data/gg2/queries.aligned.fasta",
        heldout="data/gg2/heldout.tsv",
    threads: 16
    resources:
        mem_mb=32000,
    params:
        python=sys.executable,
        script=SCRIPTS / "build_gg2.py",
        reference_dir=lambda w, input: Path(input.reference[0]).parent,
        out_dir=lambda w, output: Path(output.manifest).parent,
    log:
        "logs/gg2/build.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        """
        OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            {params.python} {params.script} --reference-dir {params.reference_dir} \
            --taxonomy {input.taxonomy} --downloads {input.downloads} \
            --out-dir {params.out_dir} --manifest {output.manifest} \
            --mem-mb {resources.mem_mb} &> {log}
        """
