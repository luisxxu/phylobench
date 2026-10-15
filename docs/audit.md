# T2 install audit

Evidence for the pins in `envs/`. Installs are proven in CI (`.github/workflows/envs.yml`), run: https://github.com/luisxxu/phylobench/actions/runs/38803340181 (all 9 install jobs, lock-guard, lint-test green). Licenses are in `LICENSES.md`.

## 1. QIIME 2 env and the `sepp-pinned` decision

- QIIME 2 2026.7 is the latest released distribution as of 2026-10-15. The distribution was renamed from "amplicon" to "qiime2" in 2026.4.
- `envs/qiime2.yaml` is a verbatim copy of `https://raw.githubusercontent.com/qiime2/distributions/4ea3dfb149410622d2d4f8219954264936111434/2026.7/qiime2/released/rachis-qiime2-linux-64-conda.yml` (fetched 2026-10-15).
- Versions in that env, from the CI log of `micromamba list`:

| Package | Version | Build |
|---|---|---|
| q2-fragment-insertion | 2026.7.0 | py312hd682af0_0 |
| sepp | 4.5.6 | py312h87e0c26_4 |
| pplacer | 1.1.alpha19 | h9ee0642_2 |
| hmmer | 3.4 | h7d74f8d_5 |
| openjdk | 25.0.2 | ha668962_0 |

- SEPP's `main.config` (via `home.path` -> `$CONDA_PREFIX/share/sepp/sepp`) resolves every tool from the env: `pplacer`, `hmmalign`, `hmmsearch`, `hmmbuild` and `seppJsonMerger.jar` all under `$CONDA_PREFIX/bin/`. SEPP uses the env's pplacer, not a bundled one.
- Rule: **"equal to our pins → no `sepp-pinned` arm"**.
- Our pins: sepp 4.5.6, pplacer 1.1.alpha23, hmmer 3.4. pplacer alpha19 ≠ alpha23, so **`sepp-pinned` is REQUIRED**.
- QIIME 2 2026.10 (dev) ships alpha23. Re-lock `envs/qiime2.yaml` when 2026.10 releases and re-evaluate this rule then (dated 2026-10-15).

## 2. pplacer 1.1.alpha23 evidence

- alpha22 -> alpha23 changes only build files and `version.ml`.
- Bioconda alpha23 is built against gsl 2.7 (bioconda-recipes PR #70019), which removes the GSL clash seen with SEPP.
- alpha19 (the QIIME 2 build) is the 2016 binary, which predates the OCaml 5 migration in alpha20. The `-j` deadlock was fixed in alpha21, so `sepp-prod` with `-x 16` is unaffected.
- Bioconda alpha22/alpha23 are built with OCaml 4.14.2.
- Observed in CI: the alpha23 build is `pplacer-1.1.alpha23-h01ede22_0`, and `pplacer --version` / `guppy --version` print only `dev`. The version proof is the conda package metadata (`micromamba list`, lockfile URL), not the binary's own output.
- The same pplacer lockfile line (URL + md5) is used by the `pplacer`, `sepp`, and `bscampp` envs; lock-guard enforces this (and the same for epa-ng in `epa-ng` and `bscampp`).

## 3. BSCAMPP config mechanism

- `bscampp/configs.py` reads `bscampp/home.path` (inside site-packages); if absent or pointing to a missing `main.config`, `bscampp/init_configs.py` regenerates `~/.bscampp/main.config`.
- `envs/bscampp.post-deploy.sh` installs bscampp by hash, writes `$CONDA_PREFIX/share/bscampp/main.config` (pplacer, epa-ng, taxit from `$CONDA_PREFIX/bin`; `hamming_distance_dir` in site-packages), and points `home.path` at that directory. Nothing is written at run time.
- CI asserts (bscampp job): `main_config_path == $CONDA_PREFIX/share/bscampp/main.config`, the three binary paths equal `$CONDA_PREFIX/bin/{pplacer,epa-ng,taxit}`, and `~/.bscampp` does not exist.
- Exception: the bundled `tools/hamming_distance/{fragment_hamming,hamming,homology}` binaries do execute. They are not deleted.
- Runtime-recompile risk: `bscampp/functions.py:ensureBinaryExecutable` treats exit code 255 (or -1) as working; anything else makes bscampp recompile into site-packages at run time. The post-deploy script fails unless all three exit 255 when run with no arguments, and they did in CI. If this check ever fails, never allow runtime recompile (T12).

## 4. APPLES-2

- Threads: `-T` (default 0 = all cores), so always pass it.
- `-s` without `-D` runs the bundled `tools/FastTree-linux`; T14 must pass `-D` with the T6 ME tree.
- `-s` also runs `TreeCluster.py` (env has `treecluster=1.0.5`).
- `apples` imports `pkg_resources`, which was removed in setuptools 82; the env pins `setuptools=81.0.0`. CI shows a `pkg_resources` deprecation warning only.
- Installed from PyPI by hash with `--no-deps --no-build-isolation`; no `wheel` pin was needed in CI.

## 5. DIPPER

- 0.1.5 has `-T/--threads` (TBB global_control; default all threads), see `dipper_cpu --help` in CI.
- `dipper_cpu --version` prints `DIPPER Version 0.1.4` for the 0.1.5 conda package; the package version is the pin.

## 6. Locks and install notes

- Lockfiles: `conda-lock==4.0.3 lock --conda <conda> -k explicit -p linux-64 -f envs/X.yaml --filename-template 'envs/X.{platform}.pin.txt'`. They contain `@EXPLICIT` and no pip lines; pip packages are installed by the `*.post-deploy.sh` scripts, which need network at env creation. Create envs on the login node (T21).
- Snakemake silently falls back from `pin.txt` to the yaml if creating from the pin fails. Wrappers must not rely on that (T8).
- glibc floor of the locked packages: see T27 for the cluster check.
- T8: Snakemake `--use-conda` calls `conda create --file`, so CI needs a conda/mamba front end for it.
- Disk: qiime2 install fit on the runner (79G free afterward); no cleanup step was needed.
