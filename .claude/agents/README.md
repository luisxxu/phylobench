# Specialist Agent Pipeline

## Mission

Six agents form the development pipeline: a Planner/SWE adversarial loop reconciles every non-trivial implementation plan before code is touched, with an explicit user-approval checkpoint, followed by QA gates and a dedicated Deployer. A read-only diagnostics agent is available independently for production troubleshooting.

## The Team

| Agent | Model | Access | Role |
|---|---|---|---|
| [`team-lead.md`](team-lead.md) | Opus | Full + dispatch | Orchestrates the pipeline; owns the final PR checklist gate |
| [`planner.md`](planner.md) | Opus | Read-only | Drafts plans; reconciles with SWE |
| [`SWE.md`](SWE.md) | Sonnet | Full | Adversarially reviews the plan, then implements it once approved |
| [`QA.md`](QA.md) | Sonnet | Full | Runs `/code-review`, full test suite, linters |
| [`deployer.md`](deployer.md) | Sonnet | Full | Ships to the target environment; owns Deploy Discipline |
| [`diagnostics.md`](diagnostics.md) | Sonnet | Read-only | Production diagnostics only, never mutates |

## Pipeline Flow

1. Team Lead triages: does this clear the trigger threshold?
2. Planner drafts a plan.
3. Planner ⇄ SWE adversarially reconcile the plan (no code written yet).
4. Team Lead surfaces the reconciled plan to the user — **hard blocking checkpoint, no implementation without explicit go-ahead.**
5. SWE implements: stubs → tests → code → `/simplify`.
6. QA runs `/code-review` → full test suite → linters.
7. Team Lead confirms the PR checklist.
8. Deployer submits HPC runs (T21/T22, and the T27 cluster check only), with a two-failed-attempt circuit breaker. Every other ticket ends at a merged PR.

`diagnostics` sits outside this flow — invoke it directly whenever you need to check the live environment's actual state.

## Trigger Threshold

Single-line/typo/comment fixes skip the pipeline entirely — just make the change. Anything touching more than one file, or any file with real invariants (schemas, migrations, config, CI, deploy manifests, secrets), goes through the full Team Lead → Planner ⇄ SWE → user go-ahead → SWE implements → QA → Deployer flow.

## Project Specifics

Placeholders are filled for phylobench: QA gates are `pytest`, `ruff`, `snakemake --lint` locally plus the CI smoke test; "deploy" means submitting SLURM runs; diagnostics is read-only `squeue`/`sacct`. Cluster hostnames are recorded in `diagnostics.md` Memory by T27. Grow each agent's `## Memory` as the project teaches you things.

## Frontmatter Convention

Every agent file in this directory starts with:

```yaml
---
name: <matches filename, case-sensitive — SWE and QA are the two uppercase exceptions>
description: <written for auto-selection — what the agent does and when to invoke it>
tools: [<restricted tool list>]
model: opus | sonnet
---
```

## Guardrails

- Do not duplicate long sprint or operational history here — link to `TICKETS.md` instead.
- Do not put secrets or customer-specific data in this directory.
- Keep `CLAUDE.md` compact — durable domain memory belongs in the relevant agent file here.
