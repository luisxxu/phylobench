# T27 cluster capability check

Barnacle (Knight Lab, SDSC), probed 2026-10-20 from an interactive job (`srun --partition=gpu -c 4 --mem=32G --gres=gpu:1 --pty bash -l`, job 1469429). No benchmark jobs were submitted.

## 1. Nodes, SLURM, storage

| Item | Value |
|---|---|
| Login node | `barnacle-login` |
| GPU node probed | `b-gpu-23` |
| CPU | 2 × Intel Xeon Silver 4110 @ 2.10 GHz, 8 cores/socket, 16 CPUs |
| RAM | 93 GiB, no swap |
| Kernel | `6.12.0-124.56.1.el10_1.x86_64` (EL10) |
| SLURM account | `knightlab` (QOS `normal`) |
| Scratch | `/ddn_scratch/$USER` (Lustre). `$SCRATCH` is not set by default. No user quota (`lfs quota` limit 0); filesystem 80% full. `/home` is NFS, 88% full — keep `data/`, `results/` and `.sif` images on scratch. |

Partitions (`sinfo -o "%P %l %G %c %m"`, one row per node class):

| Partition | Wall limit | GRES | CPUs | Memory (MB) |
|---|---|---|---|---|
| `short` (default) | 1 day | — | 64 | 514000 |
| `long` | 36 days | — | 64 | 514000 |
| `highmem` | 14 days | — | 64 | 1546528 |
| `gpu` | 14 days | `gpu:3` | 16 | 95151 |
| `gpu` | 14 days | `gpu:1` | 32 | 508000 |
| `debug` | 2 h | — | 64 | 514000 |

`qiita`, `jupyter`, `tcga`, `testvm` exist but are not for benchmark runs.

## 2. cgroup v2 `memory.peak` — yes

```
$ stat -fc %T /sys/fs/cgroup
cgroup2fs
$ cat /proc/self/cgroup
0::/system.slice/b-gpu-23_slurmstepd.scope/job_1469429/step_0/user/task_0
$ cat /sys/fs/cgroup$(cut -d: -f3 /proc/self/cgroup | head -1)/memory.peak
75055104
```

A job step can read its own task cgroup's `memory.peak`, so `mem_source=cgroup` is the primary source on Barnacle (T4, `docs/harness.md` §memory). It holds only with one rule per SLURM job (T21 profile).

## 3. GPU and Apptainer — yes

| Item | Value |
|---|---|
| GPU nodes | `b-gpu-21`, `b-gpu-22`: **RTX 3090**; `b-gpu-23`, `b-gpu-24`: **Tesla V100-PCIE-32GB** (user-reported; probe ran on `b-gpu-23`) |
| NVIDIA driver / CUDA | 580.126.20 / 13.0 |
| Apptainer | 1.5.0-1.el10_2 |
| DIPPER GPU image | `docker://swalia14/dipper@sha256:0b1fba6e…08fa2` (`envs/dipper-gpu.image.txt`) pulls and `apptainer exec --nv dipper-gpu.sif dipper --help` runs |

- Pulling the image needs ~32 GB job memory and `APPTAINER_TMPDIR`/`APPTAINER_CACHEDIR` on scratch; with less memory the extract step is OOM-killed. Compressed layers ≈ 4.3 GB.
- `--help` never touches the GPU. Whether the image runs kernels on the V100 is UNVERIFIED: CUDA 13 dropped Volta (sm_70) from the toolkit, so an image built with CUDA 13 and no sm_70 code fails on its first GPU call. T21 checks the image's CUDA version and does one real `--add` on the V100 before the pilot.
- The image's `--help` has no `-T/--threads` option (the CPU bioconda 0.1.5 build has it). The GPU arm must not pass `-T`.
- The image defaults to collapsing near-zero internal edges into multifurcations on output; `--force-binary` keeps a binary tree with `:0` branches. Relevant to the T25 RF check and GG2 backbone polytomies — decide in T21 pilot.

## 4. Snakemake SLURM executor plugin — available, not yet installed

`phylobench` env on Barnacle: `snakemake --version` → `9.27.0`; `pip show snakemake-executor-plugin-slurm` → not found. `mamba install -n phylobench --dry-run -c conda-forge -c bioconda snakemake-executor-plugin-slurm` resolves to `snakemake-executor-plugin-slurm 2.8.0 pyhdfd78af_0` + `snakemake-executor-plugin-slurm-jobstep 0.6.1 pyhdfd78af_0` (bioconda) without changing Snakemake. T21 pins these two in `envs/phylobench.yaml`.

## Open for T21

- Which card the GPU arm runs on: the RTX 3090 nodes (`b-gpu-21`/`-22`, 24 GB, design §5) or the V100 nodes (`b-gpu-23`/`-24`, 32 GB, more headroom at Bfull\* but Volta kernel support unverified). Pin it with `--nodelist` in the T21 profile; compare GPU runtimes only within one card. The 3090 nodes' CPU/RAM and the Apptainer image check there are not probed yet.
- 16 threads/job fits the `gpu:3` nodes only as a whole node (16 CPUs).
- `--force-binary` (or not) for the DIPPER GPU arm, and dropping `-T` there.
