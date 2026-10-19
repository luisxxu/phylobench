# phylobench

Benchmark of phylogenetic placement tools — pplacer (baseline), SEPP (current lab tool), pplacer-SCAMPP, EPA-ng-SCAMPP, BSCAMPP, APPLES-2, DIPPER, plus plain EPA-ng as a control — for **accuracy, runtime, and peak memory across a gradient of query-to-backbone similarity**. Goal: find a tool as accurate as the lab's **production SEPP** (production `run-sepp.sh` settings and versions on the GG2 2022.10 reference, run standalone outside QIIME 2) that is cheaper to run and suitable for eventual wide use in **Knight Lab infrastructure (Qiita / QIIME 2)**. Licenses are flagged in the report, never used to decide.

- **Design / rationale:** [`design.md`](design.md) — read the section a ticket cites before planning.
- **Backlog:** [`TICKETS.md`](TICKETS.md) — source of truth for scope and dependencies. Mirrored as GitHub issues in the `phylobench` GitHub Project.

## Principles

### KISS — Keep It Simple
- Write the smallest code that solves the problem.
- No speculative features, flags, abstractions, or fallbacks.
- Three similar lines beat a premature abstraction.
- Delete unused code; don't preserve it "just in case."

### CRUD — One Job Per Unit
- Each function does one thing: **C**reate, **R**ead, **U**pdate, or **D**elete.
- No mixed read+write side effects in a single path.
- Pure where possible; push I/O to the boundaries (Snakemake rules do I/O; `src/phylobench/` is pure and tested).
- State changes are explicit and localized.

### Style
- Code speaks for itself — no comments explaining *what*.
- Comment only non-obvious *why* (constraints, invariants, gotchas).
- Names carry meaning. If a comment is needed, try a better name first.

### Secrets
- **Never** paste auth tokens, API keys, passwords, private keys, session cookies, OAuth/JWT tokens, `.env` contents, or any post-auth credential into the chat.
- If a secret is needed, reference it by name (e.g. `SLURM_ACCOUNT`) and load it from the environment or a gitignored `.env` file.
- If a secret is accidentally shared, treat it as compromised: rotate immediately, then scrub history.
- Do not commit secrets. `.env`, `*.pem`, `*.key`, and credential files are gitignored.

## Benchmark Invariants

These are what make results comparable. Violating one invalidates a run — raise it in plan review instead.

- **Platform:** x86_64 Linux only (CI: GitHub Actions `ubuntu-latest`; production: HPC SLURM). No ARM/macOS-specific code paths; don't try to run tools natively on a Mac.
- **Pinned versions:** every tool has one env YAML in `envs/` with an exact version (see `design.md` §1.1) and a lockfile. No `latest`, no floating tags, Docker/Apptainer images pinned by digest.
- **No hidden binaries:** BSCAMPP/SCAMPP must use the pinned external `pplacer`/`epa-ng`, never their bundled copies, with an immutable per-env config written at env deploy (never `~/.bscampp`). All pplacer-based arms except `sepp-prod` share one pplacer binary. SEPP gets per-job working and temp dirs.
- **Shared inputs:** two backbones (Bfull\*, B10k\*) pruned from the GG2 2022.10 SEPP reference, shared by all arms; model numbers are the shipped ones, never re-estimated (`design.md` §1.2, §2). All aligned-input arms get the same `hmmalign` alignment.
- **Measure everything through the harness:** every tool call goes through `phylobench-measure` (T4). Peak memory = cgroup v2 `memory.peak` when readable, else psutil process-tree RSS sum; `mem_source` recorded on every row; `/usr/bin/time -v` RSS is never the primary. Phases (`ref-prep`, `align`, `place`, `convert`, `evaluate`) are measured separately.
- **Fixed threads:** 16 per job in production, set via both SLURM and each tool's flag.
- **Failures are data:** OOM/timeout/crash → `status=DNF` row. Never drop or retry silently.
- **Truth labels are honest:** truth = held-out GG2 clades/singletons. EMP ASV (and any SARS-CoV-2 MAT) positions are not ground truth (DEPP / UShER circularity) → label "concordance". GG2 truth was seeded by uDance (APPLES-2) — always caveat apples2 results.
- **Never pool similarity bins across backbone sizes.**
- **Data stays out of git:** `data/` and `results/` are gitignored; only `tests/fixtures/` is committed (<1 MB). Downloads are scripted with checksums.
- **Licenses are flagged, never decided on:** `LICENSES.md` records each arm's own license and the engine actually executed (e.g. BSCAMPP → EPA-ng AGPL-3). No license analysis; licenses never change scope, ranking, or the recommendation.

## Layout

```
workflow/Snakefile, workflow/rules/*.smk, workflow/scripts/   # orchestration (I/O)
src/phylobench/                                                # pure Python: measure, jplace, accuracy, similarity, regions, heldout, dipper_adapter, models
envs/<tool>.yaml                                               # one pinned env per tool
config/{schema,smoke,full}.yaml, config/arms/*.yaml           # run configs; one file per arm (wrappers add files only)
profiles/slurm/                                                # HPC profile (T21)
tests/, tests/fixtures/tiny/                                   # pytest + smoke fixture
docs/tools/<tool>.md, LICENSES.md                              # per-tool command/quirks, licenses
```

## Commands

```sh
mamba env create -f envs/phylobench.yaml && conda activate phylobench && pip install -e . --no-deps   # dev env
pytest                                                                   # unit tests (local + CI)
ruff check . && ruff format --check . && snakemake --lint                # lint (local + CI)
snakemake --use-conda --configfile config/smoke.yaml -c 2                # smoke test — CI (linux) only; agents on macOS check it via gh pr checks
snakemake --use-conda --profile profiles/slurm --configfile config/full.yaml   # full benchmark (HPC only, user-confirmed)
```

## Agents

A six-agent pipeline lives in `.claude/agents/` — see [`.claude/agents/README.md`](.claude/agents/README.md) for the roster, flow, and trigger threshold.

Start with `team-lead` for any change broader than a single-line fix. It dispatches Planner ⇄ SWE to reconcile a plan, blocks on your explicit go-ahead, then runs SWE → QA → Deployer. In this repo "Deployer" means **submitting benchmark runs to the HPC cluster** (T21/T22, plus the T27 cluster check); there is no service deployment.

## Workflow (One Ticket Per PR)

1. **Branch** — `git checkout -b t{N}-short-description`
2. **Plan** — Read ticket + cited `design.md` sections, outline types/functions/tests, get approval.
3. **Tests** — Write failing tests first; define expected behavior (on `tests/fixtures/tiny/`).
4. **Implement** — Minimum code to pass tests.
5. **Simplify** — `/simplify`; locally `pytest` + lint green, no warnings. Smoke test runs in CI only.
6. **PR** — `gh pr create --title "T{N}: {Title}" --body '...'` (body links the issue: `Closes #{issue}`) then `gh pr merge --auto --squash --delete-branch`.
7. **Monitor** — `gh pr checks --wait`; rebase if needed. Auto-merge fires on green.
8. **Close out** — move the ticket to **## Done** in `TICKETS.md`.

### PR Checklist (Team Lead gate)
- [ ] Ticket deliverables and acceptance criteria met; `Closes #{N}` in body
- [ ] Benchmark Invariants respected (or deviation approved in plan review)
- [ ] Tests added; local gates green; CI (incl. smoke) green
- [ ] `docs/tools/<tool>.md` updated for any tool touched; `LICENSES.md` still accurate
- [ ] No data, results, or secrets committed
- [ ] Ticket moved to **## Done**

## Formats

### Commit
```
T{N}: {verb} {what changed}
```

### PR Title
```
T{N}: {Ticket title}
```

### PR Body
Use **single quotes** around `gh pr create --body` content to avoid shell interpolation.
