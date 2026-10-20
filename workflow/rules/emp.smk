from phylobench.emp import EMP_QUERY_N

EMP_BIOM_URL = "https://ftp.microbio.me/emp/release1/otu_tables/deblur/emp_deblur_150bp.release1.biom"
# No checksum is published; pinned on first download (docs/datasets/emp.md).
EMP_BIOM_SHA256 = "b324e5098666839ab9747bd7d78e5c7a2f81cc3f249eaeea5ccd755b6c18196a"
EMP_BIOM = "data/emp/emp_deblur_150bp.release1.biom"
EMP_DATASETS = [d for d in config["datasets"] if d in EMP_QUERY_N]
EMP_BACKBONE_FILES = "tree.nwk|aligned-dna-sequences.fasta|raxml-info.txt"


rule emp_download:
    output:
        EMP_BIOM,
    params:
        fetch=SCRIPTS / "fetch.sh",
        url=EMP_BIOM_URL,
        sha256=EMP_BIOM_SHA256,
    log:
        "logs/emp/download.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.fetch} {params.url} {params.sha256} {output} &> {log}"


rule emp_query_sets:
    input:
        biom=EMP_BIOM,
        backbone_tree="data/gg2/Bfull/tree.nwk",
    output:
        "results/{dataset}/regions/V4-150/queries.fasta",
    wildcard_constraints:
        dataset="emp-10k|emp-100k",
    params:
        python=sys.executable,
        script=SCRIPTS / "emp_query_sets.py",
        n=lambda wildcards: EMP_QUERY_N[wildcards.dataset],
        seed=0,
    log:
        "results/{dataset}/regions/V4-150/emp_query_sets.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --biom {input.biom} "
        "--backbone-tree {input.backbone_tree} --n {params.n} --seed {params.seed} "
        "--out {output} &> {log}"


# EMP shares GG2's backbones, so there is one ref-prep and one HMM per backbone.
rule emp_backbone:
    input:
        "data/gg2/{backbone}/{name}",
    output:
        "data/{dataset}/{backbone}/{name}",
    wildcard_constraints:
        dataset="emp-10k|emp-100k",
        name=EMP_BACKBONE_FILES,
    log:
        "logs/emp/{dataset}/{backbone}/{name}.link.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "ln -sfr {input} {output} &> {log}"


rule emp_refprep:
    input:
        "results/ref-prep/gg2/{backbone}/{name}",
    output:
        "results/ref-prep/{dataset}/{backbone}/{name}",
    wildcard_constraints:
        dataset="emp-10k|emp-100k",
        name="me-tree.nwk|epa-ng.model",
    log:
        "logs/emp/{dataset}/{backbone}/{name}.refprep-link.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "ln -sfr {input} {output} &> {log}"


rule emp_hmm:
    input:
        "results/gg2/{backbone}/align/backbone.hmm",
    output:
        "results/{dataset}/{backbone}/align/backbone.hmm",
    wildcard_constraints:
        dataset="emp-10k|emp-100k",
    log:
        "logs/emp/{dataset}/{backbone}/hmm-link.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "ln -sfr {input} {output} &> {log}"


rule emp_manifest:
    input:
        EMP_BIOM,
        "data/gg2/MANIFEST.tsv",
        expand(
            "results/{dataset}/regions/V4-150/queries.fasta",
            dataset=EMP_DATASETS,
        ),
        expand(
            "results/{dataset}/{backbone}/V4-150/similarity.tsv",
            dataset=EMP_DATASETS,
            backbone=config["backbones"],
        ),
        expand(
            "data/{dataset}/{backbone}/{name}",
            dataset=EMP_DATASETS,
            backbone=config["backbones"],
            name=EMP_BACKBONE_FILES.split("|"),
        ),
    output:
        "data/emp/MANIFEST.tsv",
    params:
        python=sys.executable,
        script=SCRIPTS / "emp_manifest.py",
        datasets=" ".join(EMP_DATASETS),
        backbones=" ".join(config["backbones"]),
    log:
        "logs/emp/manifest.log",
    conda:
        "../../envs/phylobench.yaml"
    shell:
        "{params.python} {params.script} --datasets {params.datasets} "
        "--backbones {params.backbones} --out {output} &> {log}"
