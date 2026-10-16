import sys
from pathlib import Path


rule guppy_tog:
    input:
        jplace="{prefix}.jplace",
    output:
        tree="{prefix}.tog.tre",
        run_json="{prefix}.tog.run.json",
    log:
        "{prefix}.tog.log",
    params:
        harness=f"{sys.executable} -m phylobench.measure",
        env_lock=Path(workflow.basedir).parent / "envs" / "pplacer.linux-64.pin.txt",
    threads: config["threads"]
    conda:
        "../../envs/pplacer.yaml"
    shell:
        """
        {params.harness} --out {output.run_json} --threads {threads} \
            --tool-version 1.1.alpha23 --env-lock {params.env_lock} -- \
            guppy tog -o {output.tree} {input.jplace} > {log} 2>&1
        if [ -f {output.run_json} ] && ! grep -q '"status": "ok"' {output.run_json}; then
            [ -e {output.tree} ] || : > {output.tree}
        fi
        """
