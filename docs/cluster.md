# T27 cluster capability check

Barnacle (Knight Lab, SDSC), probed 2026-10-20 from an interactive job (`srun --partition=gpu -c 4 --mem=32G --gres=gpu:1 --pty bash -l`, job 1469429). No benchmark jobs were submitted.

## 1. Nodes, SLURM, storage

| Item | Value |
|---|---|
| Login node | `barnacle-login` |
| V100 node probed | `b-gpu-23`: 2 × Intel Xeon Silver 4110 @ 2.10 GHz, 16 CPUs, 93 GiB, no swap (`gpu:3` class) |
| RTX 3090 node probed | `b-gpu-21`: AMD EPYC 7302 16-Core, 32 CPUs, 503 GiB, no swap; RTX 3090 24576 MiB (`gpu:1` class) |
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
| GPU nodes | `b-gpu-21`, `b-gpu-22`: **RTX 3090 24 GB**; `b-gpu-23`, `b-gpu-24`: **Tesla V100-PCIE-32GB** (probed `b-gpu-21` and `b-gpu-23`) |
| NVIDIA driver / CUDA | 580.126.20 / 13.0 |
| Apptainer | 1.5.0-1.el10_2 |
| DIPPER GPU image | `docker://swalia14/dipper@sha256:0b1fba6e…08fa2` (`envs/dipper-gpu.image.txt`) pulls; `apptainer exec --nv dipper-gpu.sif dipper --help` runs on `b-gpu-21` and `b-gpu-23` |

- Pulling the image needs ~32 GB job memory and `APPTAINER_TMPDIR`/`APPTAINER_CACHEDIR` on scratch; with less memory the extract step is OOM-killed. Compressed layers ≈ 4.3 GB.
- Fixture `--add` runs on both cards (2026-10-20): `apptainer exec --nv dipper-gpu.sif dipper --add -i m -d 2 -t tests/fixtures/tiny/tree.nwk -I tiny.in.fasta -O out.nwk` (backbone + aligned queries concatenated) → `Read 110 sequences`, `Tree loaded successfully with 199 nodes`, `exit=0`, 2615-byte Newick on `b-gpu-21` (RTX 3090) and `b-gpu-23` (V100). No CUDA "no kernel image" error on the V100 despite CUDA 13 dropping Volta from the toolkit.
- The fixture is too small to show the GPU was used (`Distance Operation Time 0 ms`). The T21 pilot confirms it with `gpu_peak_mb` > 0 on a real backbone.
- The image's `--help` has no `-T/--threads` option (the CPU bioconda 0.1.5 build has it). The GPU arm must not pass `-T`.
- The image defaults to collapsing near-zero internal edges into multifurcations on output; `--force-binary` keeps a binary tree with `:0` branches. Relevant to the T25 RF check and GG2 backbone polytomies — decide in T21 pilot.

## 4. Snakemake SLURM executor plugin — available, not yet installed

`phylobench` env on Barnacle: `snakemake --version` → `9.27.0`; `pip show snakemake-executor-plugin-slurm` → not found. `mamba install -n phylobench --dry-run -c conda-forge -c bioconda snakemake-executor-plugin-slurm` resolves to `snakemake-executor-plugin-slurm 2.8.0 pyhdfd78af_0` + `snakemake-executor-plugin-slurm-jobstep 0.6.1 pyhdfd78af_0` (bioconda) without changing Snakemake. T21 pins these two in `envs/phylobench.yaml`.

## Open for T21

- Which card the GPU arm runs on: the RTX 3090 nodes (`b-gpu-21`/`-22`, 24 GB, design §5) or the V100 nodes (`b-gpu-23`/`-24`, 32 GB, more headroom at Bfull\*). Both run the image. The 3090 nodes have 32 CPUs / 503 GiB vs 16 / 93 GiB, so a 16-thread job does not need the whole node. Pin it with `--nodelist` in the T21 profile; compare GPU runtimes only within one card.
- 16 threads/job fits the V100 nodes only as a whole node (16 CPUs); the 3090 nodes have 32 CPUs.
- `--force-binary` (or not) for the DIPPER GPU arm, and dropping `-T` there.
