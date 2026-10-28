# phylobench

Benchmark of phylogenetic placement tools (pplacer, SEPP, pplacer-SCAMPP, EPA-ng-SCAMPP, BSCAMPP, APPLES-2, DIPPER, plus an EPA-ng control) for accuracy, runtime, and peak memory across a gradient of query-to-backbone sequence similarity.

- [`design.md`](design.md): design, protocol decisions, open questions
- [`TICKETS.md`](TICKETS.md): ticket backlog (mirrored as GitHub issues and a GitHub Project)
- [`CLAUDE.md`](CLAUDE.md): agent rules and benchmark invariants
- [`.claude/agents/`](.claude/agents/): agent pipeline (team-lead, planner, SWE, QA, deployer, diagnostics)

## Setup

The commands below run on x86_64 Linux. The first three commands work after T1; the smoke test works after T8.

```sh
mamba env create -f envs/phylobench.yaml && conda activate phylobench
pip install -e . --no-deps
pytest
ruff check . && ruff format --check . && snakemake --lint
snakemake --use-conda --configfile config/smoke.yaml -c 2
python workflow/scripts/plot_accuracy.py --accuracy results/accuracy.tsv --outdir results/figures && python workflow/scripts/plot_cost.py --runs results/runs.tsv --accuracy results/accuracy.tsv --outdir results/figures   # plot real tables; `snakemake plots -c 1` plots the sample tables
python workflow/scripts/report_tables.py --results results --paste results.md   # regenerate the tables in results.md from the real tables (T24); `PHYLOBENCH_RESULTS=results pytest tests/test_report.py` checks they are current
```
