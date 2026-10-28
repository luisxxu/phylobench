# phylobench results: SEPP replacement benchmark

Generated tables below are pasted by `workflow/scripts/report_tables.py` from the checksum-verified result tables of the T22 run (archive git `8fb2ab8`, `/ddn_scratch/lxxu/phylobench-archive/t22-final.tgz`, sha256 `234bbbba27ab0b3d609e6690d24a95a4eb0acc31985ce7686a7ec4016128cb37`). The tables are small aggregates of the data, not the data; the raw tables stay out of git. Every number in the prose is read from a table with its number. `tests/test_report.py` fails if this file's table blocks or the checksums below differ from what the script produces from `results/`.

Result tables (sha256, as in `results/SHA256SUMS`):

- `412956a6280818319a9067bb6275736e70557e835b1817aebc06276252b296b6` `runs.tsv`
- `4fa33d9af2af14aae241363f56393957d289eeb9fed40f0849c1a4d65ae2638a` `accuracy.tsv`
- `fefde96905f81858711de247ba0d474b0deb72405ba1f384d4835a0d033b024b` `concordance.tsv`

Regenerate the tables in this file (needs the three tables in `results/`, or point `PHYLOBENCH_RESULTS` / `--results` at them):

```sh
PYTHONPATH=src python workflow/scripts/report_tables.py --results results --paste results.md
PHYLOBENCH_RESULTS=results PYTHONPATH=src pytest tests/test_report.py
```

Figures (`results/figures/<name>.png`) are not committed; they are illustrative and every claim here rests on a table. Rebuild them with the plot command in `README.md`.

## 1. Summary

- At the production reference size (Bfull\*, about 320k leaves) no arm passes the pre-registered gate "as accurate as production SEPP" (Table 5). `sepp-pinned`, which runs the same SEPP with a newer pplacer, also fails all nine regions there, so at these sample sizes the gate cannot be passed even by a near-identical tool (Tables 3 and 5). This report calls no arm equivalent to SEPP at Bfull\*; it gives effect sizes with confidence intervals (Section 3).
- At B10k\* (10k leaves), `bscampp` and `epa-ng-scampp` pass the gate in 5 of the 6 regions with at least 200 queries and `pplacer-scampp` in 4 of 6 (Table 3); the simulated cross-check agrees (Table 6). `apples2`, `dipper-m`, `dipper-m-gpu` and `dipper-r` pass none of the gg2 regions at either size (Tables 3, 4, 5).
- At Bfull\* full-length, the three SCAMPP-family arms have about 2 edges more mean error than SEPP, with intervals that exclude 0 (Table 5).
- Cost: at B10k\* on the same CPU, `bscampp` takes 0.15 to 0.23 of `sepp-prod`'s end-to-end time (Tables 13, 15, 17, 19). Its memory is 0.20 of SEPP's at full-length and 0.72 to 0.80 in the other B10k\* cells. At Bfull\* the timing comparison with SEPP is confounded by hardware (Section 5).
- The recommendation is in Section 7. Licenses are flagged in Section 9 and are not a factor.

## 2. Methods

Design: `design.md` sections 1 to 4, 9. Run record: `docs/runbook-t22.md`.

- **Arms and pins** (`design.md` 1.1, `LICENSES.md`): `sepp-prod` (SEPP 4.5.6 with pplacer 1.1.alpha19, `-A 1000 -P 5000`, the production reference), `sepp-pinned` (same SEPP, pplacer 1.1.alpha23; a version-drift control, not a challenger), `pplacer` 1.1.alpha23, `pplacer-scampp`, `epa-ng-scampp` and `bscampp` (bscampp 1.0.8), `epa-ng` 0.3.8, `apples2` 2.0.11, `dipper-m` (CPU 0.1.5; GPU image `dipper-m-gpu`), `dipper-r` 0.1.5. `epa-ng`, `dipper-r` and `pplacer` run at B10k\* only (`pplacer` has one capped Bfull\* cell); they are controls.
- **Backbones**: B10k\* (10k leaves) and Bfull\* (the GG2 2022.10 SEPP reference minus held-out queries), shared by all arms. Model numbers are the shipped ones, never re-estimated. Bins are never pooled across backbones.
- **Alignment**: every aligned-input arm gets the same single-profile `hmmalign` alignment; SEPP aligns internally. The oracle alignment is a B10k\* diagnostic only.
- **Truth and regions**: gg2 held-out clades and singletons, nine regions (full-length, V4-150, V3-V4, V1-V2, V1-V3, V4-V5, V5-V6, V6-V8, V7-V9). The queries per region are in the `n` row of Table 4: 1000, 879, 887, 153, 144, 917, 794, 921, 46. The 10k simulated panel `sim10k` (full-length and V4-150, B10k\* only) is the non-uDance cross-check. EMP ASVs (10k and 100k) have no truth and are used for cost and concordance with `sepp-prod` only.
- **Verdict**: per backbone and region, the paired per-query difference in edge distance (tool minus `sepp-prod`), 95% bootstrap CI of the mean; pass if the upper limit is at most 0.25 edges. The tables add a descriptive reading of whether the CI is above 0 (`worse`), below 0 (`better`) or spans 0; it is not a gate. Per-bin results are descriptive.
- **Measurement**: every tool call goes through `phylobench-measure`; 16 threads; peak memory is cgroup v2 `memory.peak`; failures are `DNF` rows. The run is 304 cells and 839 `runs.tsv` rows (Table 1), about 60 node-hours; SEPP place rules at Bfull\* ran on `highmem` (Intel Xeon Gold 6242); every other place rule, and SEPP at B10k\*, ran on `long` (AMD EPYC 7302); the GPU arm ran on an RTX 3090 node (Table 2).

<!-- table:runs-overview -->
**Table 1.** Rows of runs.tsv per phase.

| phase | rows | ok | DNF | cgroup | psutil |
|---|---|---|---|---|---|
| align | 27 | 27 | 0 | 10 | 17 |
| align+place | 24 | 24 | 0 | 24 | 0 |
| align+place+convert | 24 | 24 | 0 | 24 | 0 |
| convert | 201 | 171 | 30 | 119 | 82 |
| evaluate | 304 | 304 | 0 | 304 | 0 |
| place | 256 | 256 | 0 | 256 | 0 |
| ref-prep | 3 | 3 | 0 | 3 | 0 |
<!-- /table -->

<!-- table:hardware -->
**Table 2.** Placed cells per arm, backbone and CPU.

| arm | B10k EPYC 7302 | Bfull EPYC 7302 | Bfull Xeon 6242 |
|---|---|---|---|
| sepp-prod | 13 | 0 | 11 |
| sepp-pinned | 13 | 1 | 10 |
| pplacer | 13 | 1 | 0 |
| pplacer-scampp | 13 | 11 | 0 |
| epa-ng | 13 | 0 | 0 |
| epa-ng-scampp | 13 | 11 | 0 |
| bscampp | 13 | 11 | 0 |
| apples2 | 13 | 11 | 0 |
| dipper-m | 13 | 11 | 0 |
| dipper-m-gpu | 13 | 11 | 0 |
| dipper-r | 13 | 0 | 0 |
<!-- /table -->

## 3. Verdict against production SEPP

Table 3 counts passes. Tables 4 and 5 give the verdict and effect size per region for gg2 at B10k\* and Bfull\*; Table 6 is the simulated cross-check. Illustrative figure: `results/figures/verdict.png`.

<!-- table:pass-counts -->
**Table 3.** Regions passed per arm: all regions (regions with n >= 200).

| arm | gg2 B10k | gg2 Bfull | sim10k B10k |
|---|---|---|---|
| sepp-pinned | 8/9 (5/6) | 0/9 (0/6) | 2/2 (2/2) |
| pplacer | 3/9 (3/6) | 0/1 (0/1) | 2/2 (2/2) |
| pplacer-scampp | 4/9 (4/6) | 1/9 (0/6) | 2/2 (2/2) |
| epa-ng | 4/9 (4/6) |  | 2/2 (2/2) |
| epa-ng-scampp | 5/9 (5/6) | 0/9 (0/6) | 2/2 (2/2) |
| bscampp | 5/9 (5/6) | 0/9 (0/6) | 2/2 (2/2) |
| apples2 | 0/9 (0/6) | 0/9 (0/6) | 1/2 (1/2) |
| dipper-m | 0/9 (0/6) | 0/9 (0/6) | 1/2 (1/2) |
| dipper-m-gpu | 0/9 (0/6) | 0/9 (0/6) | 1/2 (1/2) |
| dipper-r | 0/9 (0/6) |  | 0/2 (0/2) |
<!-- /table -->

<!-- table:verdict-gg2-B10k -->
**Table 4.** Verdict vs sepp-prod, gg2 at B10k\*. Cell: pre-registered verdict (pass if the upper 95% bootstrap CI of the mean paired difference vs sepp-prod is at most 0.25 edges), mean difference [95% CI] in edges (tool minus sepp-prod; negative = tool closer), then whether the CI is above 0 (`worse`), below 0 (`better`) or spans 0. `n` = queries scored.

| arm | full-length | V4-150 | V3-V4 | V1-V2 | V1-V3 | V4-V5 | V5-V6 | V6-V8 | V7-V9 |
|---|---|---|---|---|---|---|---|---|---|
| sepp-pinned | pass +0.00 [+0.00, +0.00] spans 0 | pass -0.10 [-0.22, +0.01] spans 0 | pass -0.01 [-0.02, -0.00] better | pass +0.01 [-0.06, +0.08] spans 0 | pass +0.00 [+0.00, +0.00] spans 0 | fail +0.17 [-0.00, +0.36] spans 0 | pass +0.02 [-0.00, +0.05] spans 0 | pass +0.04 [-0.10, +0.20] spans 0 | pass -0.87 [-1.39, -0.43] better |
| pplacer | fail +0.81 [+0.48, +1.18] worse | pass -0.67 [-1.03, -0.33] better | pass -0.25 [-0.47, -0.04] better | fail +1.24 [-0.10, +2.69] spans 0 | fail +1.78 [+0.46, +3.34] worse | pass -0.19 [-0.53, +0.15] spans 0 | fail -0.02 [-0.61, +0.59] spans 0 | fail +0.74 [+0.32, +1.18] worse | fail +4.17 [+1.00, +7.41] worse |
| pplacer-scampp | pass -0.27 [-0.56, +0.00] spans 0 | pass -0.59 [-1.00, -0.22] better | pass -0.49 [-0.82, -0.15] better | fail +0.50 [-0.60, +1.58] spans 0 | fail -0.01 [-0.65, +0.58] spans 0 | fail -0.04 [-0.48, +0.39] spans 0 | pass -0.41 [-0.95, +0.12] spans 0 | fail +0.08 [-0.19, +0.35] spans 0 | fail -0.48 [-3.33, +2.28] spans 0 |
| epa-ng | fail +0.27 [+0.12, +0.42] worse | pass -0.84 [-1.23, -0.48] better | pass -0.23 [-0.49, +0.02] spans 0 | fail +1.39 [-0.07, +3.00] spans 0 | fail +0.51 [-0.24, +1.26] spans 0 | pass -0.27 [-0.56, +0.00] spans 0 | pass -0.52 [-1.02, -0.04] better | fail +0.30 [-0.00, +0.62] spans 0 | fail +0.35 [-2.24, +3.09] spans 0 |
| epa-ng-scampp | pass -0.24 [-0.53, +0.03] spans 0 | pass -0.42 [-0.83, -0.02] better | pass -0.56 [-0.88, -0.25] better | fail +0.70 [-0.42, +1.79] spans 0 | fail +0.15 [-0.54, +0.79] spans 0 | pass -0.20 [-0.60, +0.19] spans 0 | pass -0.40 [-0.95, +0.13] spans 0 | fail +0.10 [-0.17, +0.38] spans 0 | fail -0.65 [-2.78, +1.54] spans 0 |
| bscampp | pass -0.18 [-0.46, +0.10] spans 0 | pass -0.43 [-0.84, -0.05] better | pass -0.55 [-0.87, -0.24] better | fail +0.53 [-0.56, +1.58] spans 0 | fail +0.14 [-0.56, +0.79] spans 0 | pass -0.23 [-0.63, +0.15] spans 0 | pass -0.46 [-1.01, +0.07] spans 0 | fail +0.14 [-0.14, +0.42] spans 0 | fail -0.20 [-3.02, +2.54] spans 0 |
| apples2 | fail +1.49 [+1.24, +1.75] worse | fail +0.71 [+0.23, +1.20] worse | fail +1.12 [+0.72, +1.53] worse | fail +5.37 [+3.80, +6.95] worse | fail +6.26 [+4.62, +7.97] worse | fail +1.95 [+1.31, +2.61] worse | fail +0.38 [-0.25, +1.00] spans 0 | fail +8.23 [+7.47, +9.00] worse | fail +6.50 [+1.98, +11.50] worse |
| dipper-m | fail +7.19 [+6.59, +7.78] worse | fail +39.55 [+38.10, +41.01] worse | fail +39.00 [+37.62, +40.45] worse | fail +37.12 [+34.12, +40.20] worse | fail +35.72 [+32.71, +38.81] worse | fail +40.97 [+39.61, +42.31] worse | fail +36.20 [+34.64, +37.73] worse | fail +41.17 [+39.87, +42.48] worse | fail +37.74 [+31.33, +44.13] worse |
| dipper-m-gpu | fail +7.21 [+6.61, +7.81] worse | fail +39.06 [+37.61, +40.55] worse | fail +38.23 [+36.84, +39.65] worse | fail +37.12 [+34.12, +40.20] worse | fail +37.51 [+34.70, +40.40] worse | fail +38.41 [+37.02, +39.82] worse | fail +38.93 [+37.32, +40.51] worse | fail +40.40 [+39.04, +41.77] worse | fail +37.11 [+30.98, +43.26] worse |
| dipper-r | fail +20.43 [+19.51, +21.37] worse | fail +19.75 [+18.54, +20.94] worse | fail +17.20 [+16.20, +18.20] worse | fail +24.13 [+21.54, +26.82] worse | fail +29.33 [+26.53, +32.26] worse | fail +19.32 [+18.23, +20.43] worse | fail +21.36 [+20.14, +22.58] worse | fail +23.71 [+22.62, +24.82] worse | fail +36.85 [+30.33, +43.33] worse |
| n | 1000 | 879 | 887 | 153 | 144 | 917 | 794 | 921 | 46 |
<!-- /table -->

<!-- table:verdict-gg2-Bfull -->
**Table 5.** Verdict vs sepp-prod, gg2 at Bfull\*. Cell: pre-registered verdict (pass if the upper 95% bootstrap CI of the mean paired difference vs sepp-prod is at most 0.25 edges), mean difference [95% CI] in edges (tool minus sepp-prod; negative = tool closer), then whether the CI is above 0 (`worse`), below 0 (`better`) or spans 0. `n` = queries scored.

| arm | full-length | V4-150 | V3-V4 | V1-V2 | V1-V3 | V4-V5 | V5-V6 | V6-V8 | V7-V9 |
|---|---|---|---|---|---|---|---|---|---|
| sepp-pinned | fail -0.05 [-0.72, +0.60] spans 0 | fail +0.53 [-1.06, +2.24] spans 0 | fail +0.11 [-0.36, +0.56] spans 0 | fail +0.48 [-0.22, +1.63] spans 0 | fail +1.35 [-0.31, +4.44] spans 0 | fail +0.59 [-1.03, +2.19] spans 0 | fail +0.65 [-0.34, +1.76] spans 0 | fail -0.06 [-1.24, +1.07] spans 0 | fail +24.63 [-12.91, +64.35] spans 0 |
| pplacer |  | unplaced |  |  |  |  |  |  |  |
| pplacer-scampp | fail +1.95 [+1.23, +2.66] worse | fail +0.52 [-1.15, +2.26] spans 0 | fail +0.48 [-1.37, +2.29] spans 0 | fail +1.45 [-2.84, +6.15] spans 0 | fail +1.99 [-1.47, +6.33] spans 0 | fail +0.64 [-0.41, +1.69] spans 0 | fail +2.57 [-0.41, +5.13] spans 0 | fail -0.18 [-3.46, +2.95] spans 0 | pass -39.33 [-85.16, -6.04] better |
| epa-ng-scampp | fail +2.22 [+1.40, +3.08] worse | fail +0.88 [-0.83, +2.69] spans 0 | fail +0.38 [-1.43, +2.21] spans 0 | fail +0.64 [-3.19, +4.65] spans 0 | fail +2.59 [-0.83, +6.94] spans 0 | fail -0.25 [-2.26, +1.60] spans 0 | fail +1.52 [-1.35, +3.92] spans 0 | fail -0.22 [-2.76, +1.80] spans 0 | fail -9.46 [-57.26, +31.65] spans 0 |
| bscampp | fail +2.26 [+1.44, +3.12] worse | fail -0.12 [-1.60, +1.38] spans 0 | fail +0.52 [-1.29, +2.32] spans 0 | fail +0.60 [-3.22, +4.61] spans 0 | fail +2.38 [-1.16, +6.78] spans 0 | fail -0.57 [-2.58, +1.30] spans 0 | fail +1.70 [-1.21, +4.09] spans 0 | fail +0.19 [-3.02, +3.22] spans 0 | fail +5.80 [-44.22, +48.89] spans 0 |
| apples2 | fail +9.36 [+8.03, +10.70] worse | fail +7.90 [+4.57, +11.75] worse | fail +3.80 [+1.88, +5.96] worse | fail +14.31 [+8.97, +19.95] worse | fail +33.74 [+17.35, +53.74] worse | fail +5.27 [+1.59, +9.42] worse | fail +38.66 [+34.10, +43.42] worse | fail +44.13 [+39.77, +48.61] worse | fail +42.33 [-21.20, +102.52] spans 0 |
| dipper-m | fail +27.65 [+24.42, +31.25] worse | fail +49.77 [+46.32, +53.67] worse | fail +55.05 [+52.12, +58.32] worse | fail +62.92 [+47.50, +81.49] worse | fail +73.41 [+56.30, +93.90] worse | fail +51.67 [+49.16, +54.51] worse | fail +52.37 [+47.88, +57.15] worse | fail +52.72 [+48.48, +57.09] worse | fail +57.87 [-5.20, +117.78] spans 0 |
| dipper-m-gpu | fail +55.20 [+51.59, +59.22] worse | fail +42.00 [+38.61, +45.82] worse | fail +47.85 [+44.96, +51.10] worse | fail +54.01 [+38.86, +72.26] worse | fail +64.32 [+47.56, +84.54] worse | fail +46.46 [+44.05, +49.24] worse | fail +47.13 [+42.75, +51.80] worse | fail +50.58 [+46.33, +54.92] worse | fail +56.15 [-6.94, +116.00] spans 0 |
| n | 1000 | 879 | 887 | 153 | 144 | 917 | 794 | 921 | 46 |
<!-- /table -->

<!-- table:verdict-sim10k-B10k -->
**Table 6.** Verdict vs sepp-prod, sim10k at B10k\*. Cell: pre-registered verdict (pass if the upper 95% bootstrap CI of the mean paired difference vs sepp-prod is at most 0.25 edges), mean difference [95% CI] in edges (tool minus sepp-prod; negative = tool closer), then whether the CI is above 0 (`worse`), below 0 (`better`) or spans 0. `n` = queries scored.

| arm | full-length | V4-150 |
|---|---|---|
| sepp-pinned | pass +0.00 [+0.00, +0.00] spans 0 | pass +0.00 [+0.00, +0.00] spans 0 |
| pplacer | pass +0.00 [+0.00, +0.00] spans 0 | pass +0.00 [+0.00, +0.01] spans 0 |
| pplacer-scampp | pass +0.00 [+0.00, +0.00] spans 0 | pass +0.00 [+0.00, +0.01] spans 0 |
| epa-ng | pass +0.03 [+0.02, +0.04] worse | pass +0.18 [+0.12, +0.24] worse |
| epa-ng-scampp | pass +0.01 [+0.01, +0.02] worse | pass +0.09 [+0.05, +0.12] worse |
| bscampp | pass +0.06 [+0.04, +0.08] worse | pass +0.14 [+0.10, +0.18] worse |
| apples2 | pass +0.07 [+0.05, +0.08] worse | fail +0.49 [+0.41, +0.58] worse |
| dipper-m | pass +0.20 [+0.17, +0.23] worse | fail +25.32 [+24.99, +25.66] worse |
| dipper-m-gpu | pass +0.20 [+0.17, +0.23] worse | fail +25.23 [+24.84, +25.61] worse |
| dipper-r | fail +2.64 [+2.32, +2.95] worse | fail +12.43 [+12.20, +12.66] worse |
| n | 1025 | 1025 |
<!-- /table -->

**B10k\*.** `bscampp` and `epa-ng-scampp` pass at full-length, V4-150, V3-V4, V4-V5 and V5-V6 and fail V6-V8 (`bscampp` +0.14 [-0.14, +0.42], `epa-ng-scampp` +0.10 [-0.17, +0.38]; Table 4). `pplacer-scampp` also fails V4-V5 (-0.04 [-0.48, +0.39]). The controls `epa-ng` and `pplacer` fail full-length (+0.27 [+0.12, +0.42] and +0.81 [+0.48, +1.18]). `sepp-pinned` passes eight of nine regions; it fails V4-V5 (+0.17 [-0.00, +0.36]) although it is SEPP with a newer pplacer. V1-V2, V1-V3 and V7-V9 (n = 153, 144, 46) are too small to decide anything: every SCAMPP-family arm and `epa-ng` fails them with an interval that spans 0 and is about 2 to 6 edges wide, and `pplacer`, `apples2` and the DIPPER arms fail them with intervals above 0 or spanning 0 (Table 4). In the simulated cross-check (Table 6) every arm except `apples2` and the DIPPER arms passes both regions; `apples2` passes full-length (+0.07) and fails V4-150 (+0.49 [+0.41, +0.58]), and `dipper-m` passes full-length (+0.20) and fails V4-150 (+25.32). `apples2` fails every gg2 region at B10k\* even though GG2 truth was seeded by uDance (APPLES-2), which would favour it (Table 4).

**Bfull\*.** No arm passes the gate in a region with at least 200 queries (Table 5). The one pass, `pplacer-scampp` at V7-V9 (n = 46, -39.33 [-85.16, -6.04]), is an artefact of a very small, very noisy region: `sepp-pinned` in the same cell is +24.63 [-12.91, +64.35]. `sepp-pinned` is the noise floor: it fails all nine regions with intervals that span 0 and are up to several edges wide (V4-150 +0.53 [-1.06, +2.24]), so a fail at Bfull\* is weak evidence about a challenger. At full-length (n = 1000), where the noise floor is -0.05 [-0.72, +0.60], `bscampp` is +2.26 [+1.44, +3.12], `epa-ng-scampp` +2.22 [+1.40, +3.08] and `pplacer-scampp` +1.95 [+1.23, +2.66]; `apples2` is +9.36, `dipper-m` +27.65 and `dipper-m-gpu` +55.20. In the other regions with at least 200 queries the SCAMPP-family intervals span 0 and are 2 to 6 edges wide, so most of those cells neither show a difference nor exclude one of about 2 edges. `pplacer` produced placements for one Bfull\* cell (V4-150) but none could be scored (`unplaced`, Section 8). The mean is carried by a tail: at full-length the SCAMPP-family medians are 1.0 to 1.5 (SEPP 1.0) and the fraction at distance 0 is 0.36 to 0.38 (SEPP 0.35; Table 7).

<!-- table:distribution-gg2 -->
**Table 7.** Edge distance per cell, gg2 full-length and V4-150 (sepp-prod included): unplaced queries, mean and median over placed queries, fraction of placed queries at distance 0.

| backbone | region | arm | queries | unplaced | mean | median | fraction at 0 |
|---|---|---|---|---|---|---|---|
| B10k | full-length | sepp-prod | 1000 | 0 | 2.70 | 1.0 | 0.35 |
| B10k | full-length | sepp-pinned | 1000 | 0 | 2.70 | 1.0 | 0.35 |
| B10k | full-length | pplacer | 1000 | 0 | 3.51 | 2.0 | 0.31 |
| B10k | full-length | pplacer-scampp | 1000 | 0 | 2.43 | 1.0 | 0.34 |
| B10k | full-length | epa-ng | 1000 | 0 | 2.97 | 1.0 | 0.26 |
| B10k | full-length | epa-ng-scampp | 1000 | 0 | 2.46 | 1.0 | 0.32 |
| B10k | full-length | bscampp | 1000 | 0 | 2.52 | 1.5 | 0.30 |
| B10k | full-length | apples2 | 1000 | 0 | 4.19 | 2.0 | 0.25 |
| B10k | full-length | dipper-m | 1000 | 0 | 9.89 | 8.0 | 0.07 |
| B10k | full-length | dipper-m-gpu | 1000 | 0 | 9.91 | 8.0 | 0.07 |
| B10k | full-length | dipper-r | 1000 | 0 | 23.13 | 21.0 | 0.02 |
| B10k | V4-150 | sepp-prod | 879 | 0 | 6.65 | 4.0 | 0.17 |
| B10k | V4-150 | sepp-pinned | 879 | 0 | 6.55 | 4.0 | 0.17 |
| B10k | V4-150 | pplacer | 879 | 0 | 5.98 | 4.0 | 0.17 |
| B10k | V4-150 | pplacer-scampp | 879 | 0 | 6.05 | 4.0 | 0.17 |
| B10k | V4-150 | epa-ng | 879 | 0 | 5.80 | 4.0 | 0.17 |
| B10k | V4-150 | epa-ng-scampp | 879 | 0 | 6.23 | 4.0 | 0.17 |
| B10k | V4-150 | bscampp | 879 | 0 | 6.21 | 4.0 | 0.17 |
| B10k | V4-150 | apples2 | 879 | 0 | 7.36 | 5.0 | 0.13 |
| B10k | V4-150 | dipper-m | 879 | 0 | 46.20 | 46.0 | 0.00 |
| B10k | V4-150 | dipper-m-gpu | 879 | 0 | 45.71 | 46.0 | 0.00 |
| B10k | V4-150 | dipper-r | 879 | 0 | 26.40 | 25.0 | 0.00 |
| Bfull | full-length | sepp-prod | 1000 | 0 | 5.26 | 1.0 | 0.35 |
| Bfull | full-length | sepp-pinned | 1000 | 0 | 5.21 | 1.0 | 0.36 |
| Bfull | full-length | pplacer-scampp | 1000 | 0 | 7.22 | 1.5 | 0.36 |
| Bfull | full-length | epa-ng-scampp | 1000 | 0 | 7.48 | 1.0 | 0.38 |
| Bfull | full-length | bscampp | 1000 | 0 | 7.52 | 1.5 | 0.38 |
| Bfull | full-length | apples2 | 1000 | 0 | 14.62 | 9.0 | 0.04 |
| Bfull | full-length | dipper-m | 1000 | 0 | 32.91 | 22.0 | 0.06 |
| Bfull | full-length | dipper-m-gpu | 1000 | 0 | 60.47 | 54.0 | 0.01 |
| Bfull | V4-150 | sepp-prod | 879 | 0 | 18.16 | 11.0 | 0.07 |
| Bfull | V4-150 | sepp-pinned | 879 | 0 | 18.68 | 11.0 | 0.07 |
| Bfull | V4-150 | pplacer | 879 | 879 | unplaced |  |  |
| Bfull | V4-150 | pplacer-scampp | 879 | 0 | 18.68 | 11.0 | 0.06 |
| Bfull | V4-150 | epa-ng-scampp | 879 | 0 | 19.04 | 11.0 | 0.06 |
| Bfull | V4-150 | bscampp | 879 | 0 | 18.04 | 11.0 | 0.06 |
| Bfull | V4-150 | apples2 | 879 | 0 | 26.06 | 15.0 | 0.05 |
| Bfull | V4-150 | dipper-m | 879 | 0 | 67.93 | 60.0 | 0.00 |
| Bfull | V4-150 | dipper-m-gpu | 879 | 0 | 60.15 | 52.0 | 0.00 |
<!-- /table -->

## 4. Per-bin descriptive results

Per-bin results are descriptive only (`design.md` section 3). Tables 9 to 12 cover gg2 full-length and V4-150 at B10k\* and Bfull\*; bins with fewer than 100 queries are marked. Table 8 is the oracle-alignment diagnostic. Illustrative figures: `results/figures/accuracy_edge_dist.png`, `results/figures/oracle_gap.png`.

- At Bfull\* full-length the excess error of the SCAMPP family sits in the 75-90% bins: 85-90% (n = 177) 12.49 to 12.64 against 4.54 for SEPP; the 75-80% and 80-85% bins have 26 and 38 queries (Table 10). From 90% up the means are close to SEPP's, except the 99-99.9% bin (n = 103, 11.98 to 15.18 against 11.42).
- At B10k\* the V4-150 means are close to each other for SEPP and the likelihood arms in the bins with at least 100 queries (Table 11); at V4-150 the fraction of placements at distance 0 is 0.17 for SEPP and every likelihood arm (Table 7), which reflects fragments that fit several edges equally.
- The oracle alignment is closer than `hmmalign` for the likelihood arms by 0.5 to 1.3 edges in V5-V6 and V6-V8 (for example `bscampp` -1.14 and -0.74; Table 8): this is consistent with part of the V6-V8 fail above coming from the single-profile alignment rather than the placement algorithm (an inference; the oracle alignment was not run in the verdict). The DIPPER arms do not follow this pattern (for example `dipper-m-gpu` +2.45 [+1.85, +3.05] at V3-V4).

<!-- table:oracle-gg2-B10k -->
**Table 8.** Oracle minus hmmalign alignment, same arm, gg2 at B10k\*: mean edge-distance difference [95% CI]; negative = oracle alignment is closer.

| arm | full-length | V4-150 | V3-V4 | V1-V2 | V1-V3 | V4-V5 | V5-V6 | V6-V8 | V7-V9 |
|---|---|---|---|---|---|---|---|---|---|
| pplacer | -0.12 [-0.20, -0.04] | -0.03 [-0.14, +0.09] | -0.11 [-0.20, -0.03] | +0.00 [-0.52, +0.52] | -0.15 [-0.51, +0.21] | -0.26 [-0.48, -0.06] | -1.10 [-1.42, -0.78] | -0.54 [-0.75, -0.33] | -0.48 [-1.50, +0.24] |
| pplacer-scampp | -0.08 [-0.19, +0.01] | -0.13 [-0.29, +0.03] | -0.12 [-0.25, +0.00] | -0.35 [-1.08, +0.40] | -0.06 [-0.54, +0.48] | -0.19 [-0.36, -0.03] | -1.16 [-1.52, -0.81] | -0.67 [-0.94, -0.42] | -1.24 [-2.72, -0.15] |
| epa-ng | -0.09 [-0.21, +0.02] | -0.05 [-0.25, +0.15] | -0.15 [-0.33, +0.03] | -1.52 [-2.99, -0.27] | -0.44 [-1.06, +0.20] | -0.17 [-0.32, -0.03] | -1.26 [-1.57, -0.94] | -0.76 [-1.05, -0.50] | -2.02 [-5.54, +1.04] |
| epa-ng-scampp | -0.12 [-0.23, -0.03] | -0.05 [-0.23, +0.13] | -0.12 [-0.25, +0.00] | -0.30 [-1.03, +0.42] | -0.18 [-0.74, +0.40] | -0.08 [-0.23, +0.07] | -1.19 [-1.51, -0.87] | -0.76 [-0.99, -0.54] | -1.52 [-3.67, +0.09] |
| bscampp | -0.12 [-0.23, -0.02] | +0.20 [-0.04, +0.47] | -0.16 [-0.32, -0.02] | -0.14 [-0.81, +0.56] | -0.16 [-0.71, +0.42] | -0.01 [-0.19, +0.18] | -1.14 [-1.47, -0.81] | -0.74 [-0.97, -0.52] | -0.57 [-3.70, +2.41] |
| apples2 | -0.31 [-0.49, -0.13] | +0.08 [-0.03, +0.20] | -0.15 [-0.29, -0.03] | -0.39 [-1.39, +0.59] | +0.29 [+0.06, +0.57] | -0.17 [-0.62, +0.25] | +0.38 [-0.13, +0.91] | -1.81 [-2.23, -1.41] | -0.65 [-1.74, +0.09] |
| dipper-m | -1.08 [-1.53, -0.64] | +0.37 [-0.09, +0.83] | -0.38 [-0.84, +0.08] | -3.48 [-5.42, -1.39] | -1.72 [-2.85, -0.53] | -3.13 [-3.71, -2.55] | +2.24 [+1.60, +2.89] | -2.02 [-2.50, -1.54] | -0.70 [-2.28, +0.91] |
| dipper-m-gpu | -1.03 [-1.49, -0.57] | -1.17 [-1.62, -0.73] | +2.45 [+1.85, +3.05] | -3.56 [-5.64, -1.40] | -2.21 [-3.81, -0.53] | -0.58 [-1.06, -0.12] | -3.64 [-4.29, -2.99] | -0.59 [-1.04, -0.14] | -0.07 [-1.24, +1.20] |
<!-- /table -->

<!-- table:bins-gg2-B10k-full-length -->
**Table 9.** Per-bin mean edge distance, gg2 full-length at B10k\*. Mean edge distance per similarity bin; `†` = fewer than 100 queries; `n` = queries in the bin. Descriptive only, never pooled across backbones.

| arm | <75 | 75-80 | 80-85 | 85-90 | 90-95 | 95-97 | 97-99 | 99-99.9 |
|---|---|---|---|---|---|---|---|---|
| sepp-prod | 2.55† | 5.22† | 6.89 | 3.36 | 2.52 | 0.83 | 1.10 | 2.30† |
| sepp-pinned | 2.55† | 5.22† | 6.89 | 3.36 | 2.52 | 0.83 | 1.10 | 2.30† |
| pplacer | 2.52† | 4.94† | 7.14 | 2.92 | 2.52 | 1.16 | 3.70 | 7.57† |
| pplacer-scampp | 3.00† | 5.44† | 4.18 | 2.82 | 2.47 | 1.02 | 1.32 | 2.98† |
| epa-ng | 2.62† | 5.89† | 7.33 | 2.98 | 2.70 | 1.21 | 1.67 | 3.15† |
| epa-ng-scampp | 3.17† | 5.42† | 4.23 | 2.99 | 2.34 | 1.06 | 1.35 | 3.19† |
| bscampp | 3.17† | 6.11† | 4.61 | 3.05 | 2.34 | 1.05 | 1.36 | 2.98† |
| apples2 | 1.97† | 7.28† | 9.03 | 3.47 | 4.06 | 2.29 | 3.24 | 4.68† |
| dipper-m | 10.66† | 11.67† | 14.83 | 6.21 | 11.02 | 8.60 | 8.84 | 11.11† |
| dipper-m-gpu | 10.93† | 11.67† | 14.83 | 6.28 | 11.02 | 8.60 | 8.86 | 11.13† |
| dipper-r | 19.24† | 14.11† | 25.51 | 20.05 | 22.51 | 23.81 | 23.72 | 37.13† |
| n | 29 | 36 | 106 | 165 | 281 | 112 | 224 | 47 |
<!-- /table -->

<!-- table:bins-gg2-Bfull-full-length -->
**Table 10.** Per-bin mean edge distance, gg2 full-length at Bfull\*. Mean edge distance per similarity bin; `†` = fewer than 100 queries; `n` = queries in the bin. Descriptive only, never pooled across backbones.

| arm | <75 | 75-80 | 80-85 | 85-90 | 90-95 | 95-97 | 97-99 | 99-99.9 |
|---|---|---|---|---|---|---|---|---|
| sepp-prod | 0.00† | 8.23† | 2.61† | 4.54 | 6.31 | 3.44† | 3.74 | 11.42 |
| sepp-pinned | 0.00† | 8.23† | 2.61† | 4.54 | 6.31 | 3.52† | 3.71 | 10.94 |
| pplacer-scampp | 0.00† | 14.15† | 8.61† | 12.49 | 6.16 | 4.12† | 4.01 | 11.98 |
| epa-ng-scampp | 0.00† | 13.88† | 8.61† | 12.64 | 6.06 | 4.15† | 3.91 | 14.89 |
| bscampp | 0.00† | 14.12† | 8.66† | 12.61 | 6.08 | 4.16† | 3.91 | 15.18 |
| apples2 | 4.10† | 26.19† | 17.71† | 10.32 | 7.29 | 8.91† | 19.93 | 25.46 |
| dipper-m | 35.40† | 22.38† | 35.95† | 28.27 | 24.75 | 28.09† | 34.02 | 61.38 |
| dipper-m-gpu | 47.85† | 21.88† | 33.97† | 31.41 | 48.69 | 56.79† | 79.65 | 104.22 |
| n | 20 | 26 | 38 | 177 | 240 | 82 | 314 | 103 |
<!-- /table -->

<!-- table:bins-gg2-B10k-V4-150 -->
**Table 11.** Per-bin mean edge distance, gg2 V4-150 at B10k\*. Mean edge distance per similarity bin; `†` = fewer than 100 queries; `n` = queries in the bin. Descriptive only, never pooled across backbones.

| arm | <75 | 75-80 | 80-85 | 85-90 | 90-95 | 95-97 | 97-99 | 99-99.9 | >=99.9 |
|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | 11.00† | 16.51† | 13.46† | 9.25† | 5.84 | 4.76 | 2.90 | 3.50† | 5.95 |
| sepp-pinned | 11.00† | 16.51† | 13.46† | 9.25† | 5.84 | 4.76 | 2.90 | 3.42† | 5.54 |
| pplacer | 16.00† | 17.04† | 9.59† | 6.44† | 5.55 | 4.60 | 2.34 | 3.01† | 6.02 |
| pplacer-scampp | 17.67† | 20.31† | 8.95† | 8.34† | 5.31 | 4.04 | 2.33 | 3.10† | 5.37 |
| epa-ng | 15.67† | 14.75† | 9.36† | 6.57† | 5.62 | 4.66 | 2.69 | 2.92† | 5.63 |
| epa-ng-scampp | 15.67† | 20.33† | 10.02† | 8.09† | 5.33 | 4.66 | 2.59 | 3.21† | 5.36 |
| bscampp | 15.67† | 20.41† | 9.47† | 8.01† | 5.33 | 4.65 | 2.59 | 3.08† | 5.51 |
| apples2 | 12.00† | 17.86† | 15.17† | 11.46† | 6.15 | 4.10 | 3.56 | 4.31† | 6.76 |
| dipper-m | 16.33† | 26.63† | 26.73† | 34.01† | 42.96 | 47.15 | 47.81 | 55.14† | 60.34 |
| dipper-m-gpu | 18.67† | 20.12† | 28.17† | 30.78† | 42.89 | 47.02 | 48.64 | 54.04† | 61.00 |
| dipper-r | 7.67† | 11.31† | 10.71† | 16.07† | 20.88 | 25.55 | 33.09 | 33.17† | 38.59 |
| n | 3 | 51 | 59 | 91 | 156 | 137 | 101 | 80 | 201 |
<!-- /table -->

<!-- table:bins-gg2-Bfull-V4-150 -->
**Table 12.** Per-bin mean edge distance, gg2 V4-150 at Bfull\*. Mean edge distance per similarity bin; `†` = fewer than 100 queries; `n` = queries in the bin. Descriptive only, never pooled across backbones.

| arm | 75-80 | 80-85 | 85-90 | 90-95 | 95-97 | 97-99 | 99-99.9 | >=99.9 |
|---|---|---|---|---|---|---|---|---|
| sepp-prod | 13.05† | 16.72† | 14.88† | 15.76† | 15.67† | 10.74 | 18.64† | 22.60 |
| sepp-pinned | 13.05† | 16.72† | 14.88† | 15.76† | 15.66† | 10.62 | 16.00† | 24.39 |
| pplacer | unplaced | unplaced | unplaced | unplaced | unplaced | unplaced | unplaced | unplaced |
| pplacer-scampp | 29.37† | 24.11† | 13.56† | 14.86† | 13.18† | 13.06 | 14.86† | 22.72 |
| epa-ng-scampp | 30.00† | 24.57† | 14.18† | 14.58† | 13.12† | 12.43 | 15.88† | 23.44 |
| bscampp | 30.00† | 24.97† | 13.46† | 14.55† | 13.02† | 12.42 | 15.58† | 21.24 |
| apples2 | 37.16† | 30.39† | 14.86† | 15.74† | 20.18† | 17.38 | 22.19† | 33.74 |
| dipper-m | 28.95† | 27.62† | 37.30† | 42.50† | 53.70† | 54.56 | 70.34† | 93.87 |
| dipper-m-gpu | 16.79† | 19.43† | 30.33† | 35.42† | 45.22† | 44.82 | 61.22† | 87.18 |
| n | 19 | 61 | 57 | 84 | 83 | 125 | 74 | 376 |
<!-- /table -->

## 5. Cost

Tables 13 to 20 give placement and end-to-end time, peak memory and the ratio to `sepp-prod` in the same cell. Illustrative figures: `results/figures/cost_vs_backbone.png`, `results/figures/cost_vs_query_count.png`, `results/figures/end_to_end.png`, `results/figures/gpu_memory.png`.

- **Same hardware (B10k\*).** All B10k\* cells ran on the EPYC nodes (Table 2), so the ratios are like for like. `bscampp` takes 0.15 (full-length), 0.16 (V4-150), 0.23 (EMP 10k) and 0.21 (EMP 100k) of `sepp-prod`'s end-to-end time (Tables 13, 15, 17, 19) and 0.20, 0.80, 0.79 and 0.72 of its memory. `epa-ng-scampp` is 0.92, 1.56, 2.60 and 0.77 in time and 0.37, 1.55, 2.77 and 3.01 in memory. `pplacer-scampp` is slower than SEPP in every B10k\* cell (1.50 to 4.84) and uses 0.40 to 2.33 of its memory. The controls: `epa-ng` 0.13 to 0.21 in time, 0.15 to 0.54 in memory; `pplacer` 0.46 to 1.33 in time, 0.14 to 0.32 in memory.
- **Bfull\* wall time is confounded.** SEPP at Bfull\* ran on Xeon 6242 nodes and every challenger on EPYC 7302 nodes (Table 2; one `sepp-pinned` cell, Table 16, ran on EPYC). No hardware-normalised claim is made. For `bscampp` the end-to-end ratio to `sepp-prod` is 0.15, 0.20, 0.73 and 0.26 across the four Bfull\* cost cells (Tables 14, 16, 18, 20; 645.6 s against 4243.7 s at full-length); `epa-ng-scampp` 0.20, 0.30, 1.00 and 0.87; `pplacer-scampp` 0.37, 0.37, 1.27 and 1.11. Only B10k\* shows the same-hardware ratios above.
- **Memory** is measured by cgroup `memory.peak` on every placement row of the cost tables and is comparable across machines. At Bfull\* full-length SEPP peaks at 277 GB and `bscampp` at 22.8 GB (Table 14); in the 100k EMP cell `epa-ng-scampp` (2.55) and `pplacer-scampp` (2.44) use more memory than SEPP while `bscampp` uses 0.33 (Table 20).
- **GPU, beside the CPU ranking.** `dipper-m-gpu` is among the cheapest placers (end-to-end 0.06 to 0.22 of `sepp-prod`; GPU memory 336 to 900 MB) but fails the accuracy gate everywhere, so it does not compete in the ranking (Tables 13 to 20).
- **Unfinished end-to-end cells.** `apples2` cannot be ranked end to end: its `guppy tog` conversion fails in every Bfull\* cell and in the EMP cells (marked `DNF (convert)`), although its placements are scored. `pplacer` at Bfull\* V4-150 placed in 9233 s with 73.8 GB (Table 16) and its conversion failed.

<!-- table:cost-gg2-B10k-full-length -->
**Table 13.** Cost, gg2 full-length at B10k\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 1000 | 337.6 | 337.6 | 1.00 | 102337 | 1.00 |  | cgroup | EPYC 7302 | ok |
| sepp-pinned | raw | 1000 | 286.9 | 287.5 | 0.85 | 103229 | 1.01 |  | cgroup | EPYC 7302 | ok |
| pplacer | hmmalign | 1000 | 248.3 | 269.2 | 0.80 | 14795 | 0.14 |  | cgroup | EPYC 7302 | ok |
| pplacer-scampp | hmmalign | 1000 | 934.0 | 954.8 | 2.83 | 41414 | 0.40 |  | cgroup | EPYC 7302 | ok |
| epa-ng | hmmalign | 1000 | 23.1 | 44.1 | 0.13 | 15811 | 0.15 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 1000 | 288.8 | 309.6 | 0.92 | 37729 | 0.37 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 1000 | 30.7 | 51.6 | 0.15 | 20943 | 0.20 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 1000 | 26.1 | 47.1 | 0.14 | 737 | 0.01 |  | cgroup | EPYC 7302 | ok |
| dipper-m | hmmalign | 1000 | 4.4 | 24.6 | 0.07 | 194 | 0.00 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 1000 | 4.8 | 25.0 | 0.07 | 585 | 0.01 | 336 | cgroup | EPYC 7302 | ok |
| dipper-r | raw | 1000 | 6.5 | 6.5 | 0.02 | 340 | 0.00 |  | cgroup | EPYC 7302 | ok |
<!-- /table -->

<!-- table:cost-gg2-Bfull-full-length -->
**Table 14.** Cost, gg2 full-length at Bfull\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 1000 | 4243.7 | 4243.7 | 1.00 | 277165 | 1.00 |  | cgroup | Xeon 6242 | ok |
| sepp-pinned | raw | 1000 | 3939.7 | 3982.1 | 0.94 | 262963 | 0.95 |  | cgroup | Xeon 6242 | ok |
| pplacer-scampp | hmmalign | 1000 | 1507.9 | 1568.1 | 0.37 | 42929 | 0.15 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 1000 | 812.3 | 869.1 | 0.20 | 40107 | 0.14 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 1000 | 588.4 | 645.6 | 0.15 | 22757 | 0.08 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 1000 | 269.6 | 319.2 |  | 15459 | 0.06 |  | cgroup | EPYC 7302 | DNF (convert) |
| dipper-m | hmmalign | 1000 | 127.6 | 147.3 | 0.03 | 2771 | 0.01 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 1000 | 219.0 | 238.7 | 0.06 | 2906 | 0.01 | 780 | cgroup | EPYC 7302 | ok |
<!-- /table -->

<!-- table:cost-gg2-B10k-V4-150 -->
**Table 15.** Cost, gg2 V4-150 at B10k\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 879 | 62.5 | 62.5 | 1.00 | 11566 | 1.00 |  | cgroup | EPYC 7302 | ok |
| sepp-pinned | raw | 879 | 58.8 | 59.3 | 0.95 | 11551 | 1.00 |  | cgroup | EPYC 7302 | ok |
| pplacer | hmmalign | 879 | 25.9 | 28.9 | 0.46 | 2482 | 0.21 |  | cgroup | EPYC 7302 | ok |
| pplacer-scampp | hmmalign | 879 | 175.6 | 177.7 | 2.84 | 7508 | 0.65 |  | cgroup | EPYC 7302 | ok |
| epa-ng | hmmalign | 879 | 10.1 | 12.2 | 0.19 | 6198 | 0.54 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 879 | 95.3 | 97.4 | 1.56 | 17951 | 1.55 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 879 | 8.1 | 10.2 | 0.16 | 9207 | 0.80 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 879 | 4.8 | 6.7 |  | 557 | 0.05 |  | cgroup | EPYC 7302 | DNF (convert) |
| dipper-m | hmmalign | 879 | 3.1 | 4.7 | 0.07 | 192 | 0.02 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 879 | 4.0 | 5.5 | 0.09 | 586 | 0.05 | 336 | cgroup | EPYC 7302 | ok |
| dipper-r | raw | 879 | 4.2 | 4.2 | 0.07 | 340 | 0.03 |  | cgroup | EPYC 7302 | ok |
<!-- /table -->

<!-- table:cost-gg2-Bfull-V4-150 -->
**Table 16.** Cost, gg2 V4-150 at Bfull\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 879 | 1263.0 | 1263.0 | 1.00 | 36017 | 1.00 |  | cgroup | Xeon 6242 | ok |
| sepp-pinned | raw | 879 | 1054.4 | 1091.5 | 0.86 | 33819 | 0.94 |  | cgroup | EPYC 7302 | ok |
| pplacer | hmmalign | 879 | 9233.1 | 9235.2 |  | 73822 | 2.05 |  | cgroup | EPYC 7302 | DNF (convert) |
| pplacer-scampp | hmmalign | 879 | 427.8 | 465.9 | 0.37 | 10675 | 0.30 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 879 | 337.8 | 379.4 | 0.30 | 19759 | 0.55 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 879 | 216.7 | 256.1 | 0.20 | 11531 | 0.32 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 879 | 95.5 | 122.5 |  | 6650 | 0.18 |  | cgroup | EPYC 7302 | DNF (convert) |
| dipper-m | hmmalign | 879 | 101.2 | 102.8 | 0.08 | 2794 | 0.08 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 879 | 215.7 | 217.3 | 0.17 | 2915 | 0.08 | 780 | cgroup | EPYC 7302 | ok |
<!-- /table -->

<!-- table:cost-emp-10k-B10k-V4-150 -->
**Table 17.** Cost, emp-10k V4-150 at B10k\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 10000 | 242.6 | 242.6 | 1.00 | 12212 | 1.00 |  | cgroup | EPYC 7302 | ok |
| sepp-pinned | raw | 10000 | 207.1 | 208.3 | 0.86 | 12133 | 0.99 |  | cgroup | EPYC 7302 | ok |
| pplacer | hmmalign | 10000 | 239.1 | 257.2 | 1.06 | 3546 | 0.29 |  | cgroup | EPYC 7302 | ok |
| pplacer-scampp | hmmalign | 10000 | 1153.3 | 1173.4 | 4.84 | 25535 | 2.09 |  | cgroup | EPYC 7302 | ok |
| epa-ng | hmmalign | 10000 | 22.8 | 41.6 | 0.17 | 6290 | 0.52 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 10000 | 613.6 | 631.7 | 2.60 | 33811 | 2.77 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 10000 | 37.4 | 55.3 | 0.23 | 9671 | 0.79 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 10000 | 42.1 | 59.9 |  | 885 | 0.07 |  | cgroup | EPYC 7302 | DNF (convert) |
| dipper-m | hmmalign | 10000 | 44.3 | 61.0 | 0.25 | 252 | 0.02 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 10000 | 11.0 | 27.7 | 0.11 | 631 | 0.05 | 348 | cgroup | EPYC 7302 | ok |
| dipper-r | raw | 10000 | 69.5 | 69.5 | 0.29 | 481 | 0.04 |  | cgroup | EPYC 7302 | ok |
<!-- /table -->

<!-- table:cost-emp-10k-Bfull-V4-150 -->
**Table 18.** Cost, emp-10k V4-150 at Bfull\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 10000 | 3865.7 | 3865.7 | 1.00 | 45103 | 1.00 |  | cgroup | Xeon 6242 | ok |
| sepp-pinned | raw | 10000 | 3653.3 | 3691.8 | 0.96 | 45565 | 1.01 |  | cgroup | Xeon 6242 | ok |
| pplacer-scampp | hmmalign | 10000 | 4869.5 | 4922.8 | 1.27 | 45135 | 1.00 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 10000 | 3816.0 | 3872.8 | 1.00 | 53929 | 1.20 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 10000 | 2750.5 | 2807.6 | 0.73 | 13437 | 0.30 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 10000 | 473.6 | 515.1 |  | 16802 | 0.37 |  | cgroup | EPYC 7302 | DNF (convert) |
| dipper-m | hmmalign | 10000 | 1013.3 | 1030.0 | 0.27 | 2799 | 0.06 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 10000 | 318.8 | 335.5 | 0.09 | 2582 | 0.06 | 786 | cgroup | EPYC 7302 | ok |
<!-- /table -->

<!-- table:cost-emp-100k-B10k-V4-150 -->
**Table 19.** Cost, emp-100k V4-150 at B10k\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 100000 | 1919.8 | 1919.8 | 1.00 | 16184 | 1.00 |  | cgroup | EPYC 7302 | ok |
| sepp-pinned | raw | 100000 | 1503.2 | 1526.5 | 0.80 | 15367 | 0.95 |  | cgroup | EPYC 7302 | ok |
| pplacer | hmmalign | 100000 | 2357.9 | 2557.1 | 1.33 | 5232 | 0.32 |  | cgroup | EPYC 7302 | ok |
| pplacer-scampp | hmmalign | 100000 | 2696.5 | 2884.8 | 1.50 | 37771 | 2.33 |  | cgroup | EPYC 7302 | ok |
| epa-ng | hmmalign | 100000 | 164.4 | 410.7 | 0.21 | 6541 | 0.40 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 100000 | 1289.7 | 1476.7 | 0.77 | 48762 | 3.01 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 100000 | 219.1 | 410.2 | 0.21 | 11706 | 0.72 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 100000 | 335.4 | 512.2 |  | 1633 | 0.10 |  | cgroup | EPYC 7302 | DNF (convert) |
| dipper-m | hmmalign | 100000 | 1667.3 | 1836.8 | 0.96 | 815 | 0.05 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 100000 | 250.0 | 419.5 | 0.22 | 776 | 0.05 | 474 | cgroup | EPYC 7302 | ok |
| dipper-r | raw | 100000 | 3295.5 | 3295.5 | 1.72 | 1916 | 0.12 |  | cgroup | EPYC 7302 | ok |
<!-- /table -->

<!-- table:cost-emp-100k-Bfull-V4-150 -->
**Table 20.** Cost, emp-100k V4-150 at Bfull\*. `place s` = placement wall time; `e2e s` = align + place + convert; ratios are to sepp-prod in the same cell and blank where the cell did not finish. Peak memory is the place phase. `gpu peak MB` is GPU memory (dipper-m-gpu only).

| arm | aln | queries | place s | e2e s | e2e vs sepp-prod | peak mem MB | mem vs sepp-prod | gpu peak MB | mem source | cpu | e2e status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sepp-prod | raw | 100000 | 22904.0 | 22904.0 | 1.00 | 47715 | 1.00 |  | cgroup | Xeon 6242 | ok |
| sepp-pinned | raw | 100000 | 22826.2 | 22880.2 | 1.00 | 48481 | 1.02 |  | cgroup | Xeon 6242 | ok |
| pplacer-scampp | hmmalign | 100000 | 25218.2 | 25438.1 | 1.11 | 116419 | 2.44 |  | cgroup | EPYC 7302 | ok |
| epa-ng-scampp | hmmalign | 100000 | 19625.6 | 19841.4 | 0.87 | 121589 | 2.55 |  | cgroup | EPYC 7302 | ok |
| bscampp | hmmalign | 100000 | 5727.1 | 5948.0 | 0.26 | 15936 | 0.33 |  | cgroup | EPYC 7302 | ok |
| apples2 | hmmalign | 100000 | 4362.2 | 4557.7 |  | 18672 | 0.39 |  | cgroup | EPYC 7302 | DNF (convert) |
| dipper-m | hmmalign | 100000 | 11335.4 | 11505.4 | 0.50 | 3072 | 0.06 |  | cgroup | EPYC 7302 | ok |
| dipper-m-gpu | hmmalign | 100000 | 1519.8 | 1689.8 | 0.07 | 3143 | 0.07 | 900 | cgroup | EPYC 7302 | ok |
<!-- /table -->

## 6. EMP concordance

Concordance is the edge distance between an arm's placement of an EMP ASV and `sepp-prod`'s. The EMP ASVs have no truth, and GG2's own placement of them came from DEPP, so this is agreement with SEPP and never accuracy. No ASV is unplaced in any arm (Table 21). The `sepp-pinned` rows are the noise floor: mean 0.01 at B10k\* and 0.32 to 0.34 at Bfull\* (fraction at 0 of 1.00 and 0.98). At B10k\* (EMP 10k) the other arms are at mean 3.88 (`pplacer`), 4.54 (`pplacer-scampp`), 4.51 (`epa-ng`), 4.82 (`epa-ng-scampp`), 4.83 (`bscampp`), 7.77 (`apples2`) and 38.25 (`dipper-m`); at Bfull\* 8.69, 9.30, 9.30, 18.19 and 78.08 for `pplacer-scampp`, `epa-ng-scampp`, `bscampp`, `apples2` and `dipper-m`. The 100k rows are within 0.1 of the 10k rows for the SCAMPP family. These rows say that the SCAMPP-family arms agree with SEPP less at Bfull\* than at B10k\*, in line with Table 5; they say nothing about which is right.

<!-- table:concordance -->
**Table 21.** EMP concordance with sepp-prod (edge distance between each arm's placement and sepp-prod's; not accuracy). Mean and median over placed ASVs.

| dataset | backbone | arm | queries | unplaced | mean | median | fraction at 0 |
|---|---|---|---|---|---|---|---|
| emp-100k | B10k | sepp-prod | 100000 | 0 | 0.00 | 0.0 | 1.00 |
| emp-100k | B10k | sepp-pinned | 100000 | 0 | 0.01 | 0.0 | 1.00 |
| emp-100k | B10k | pplacer | 100000 | 0 | 3.85 | 0.0 | 0.75 |
| emp-100k | B10k | pplacer-scampp | 100000 | 0 | 4.52 | 0.0 | 0.73 |
| emp-100k | B10k | epa-ng | 100000 | 0 | 4.53 | 0.0 | 0.61 |
| emp-100k | B10k | epa-ng-scampp | 100000 | 0 | 4.81 | 0.0 | 0.65 |
| emp-100k | B10k | bscampp | 100000 | 0 | 4.84 | 0.0 | 0.64 |
| emp-100k | B10k | apples2 | 100000 | 0 | 7.66 | 4.0 | 0.26 |
| emp-100k | B10k | dipper-m | 100000 | 0 | 37.90 | 37.0 | 0.00 |
| emp-100k | B10k | dipper-m-gpu | 100000 | 0 | 37.77 | 37.0 | 0.00 |
| emp-100k | B10k | dipper-r | 100000 | 0 | 23.49 | 23.0 | 0.04 |
| emp-100k | Bfull | sepp-prod | 100000 | 0 | 0.00 | 0.0 | 1.00 |
| emp-100k | Bfull | sepp-pinned | 100000 | 0 | 0.32 | 0.0 | 0.98 |
| emp-100k | Bfull | pplacer-scampp | 100000 | 0 | 8.63 | 0.0 | 0.68 |
| emp-100k | Bfull | epa-ng-scampp | 100000 | 0 | 9.28 | 0.0 | 0.56 |
| emp-100k | Bfull | bscampp | 100000 | 0 | 9.31 | 0.0 | 0.55 |
| emp-100k | Bfull | apples2 | 100000 | 0 | 18.38 | 12.0 | 0.16 |
| emp-100k | Bfull | dipper-m | 100000 | 0 | 89.94 | 76.0 | 0.00 |
| emp-100k | Bfull | dipper-m-gpu | 100000 | 0 | 83.93 | 82.0 | 0.00 |
| emp-10k | B10k | sepp-prod | 10000 | 0 | 0.00 | 0.0 | 1.00 |
| emp-10k | B10k | sepp-pinned | 10000 | 0 | 0.01 | 0.0 | 1.00 |
| emp-10k | B10k | pplacer | 10000 | 0 | 3.88 | 0.0 | 0.75 |
| emp-10k | B10k | pplacer-scampp | 10000 | 0 | 4.54 | 0.0 | 0.73 |
| emp-10k | B10k | epa-ng | 10000 | 0 | 4.51 | 0.0 | 0.61 |
| emp-10k | B10k | epa-ng-scampp | 10000 | 0 | 4.82 | 0.0 | 0.65 |
| emp-10k | B10k | bscampp | 10000 | 0 | 4.83 | 0.0 | 0.64 |
| emp-10k | B10k | apples2 | 10000 | 0 | 7.77 | 4.0 | 0.26 |
| emp-10k | B10k | dipper-m | 10000 | 0 | 38.25 | 37.0 | 0.00 |
| emp-10k | B10k | dipper-m-gpu | 10000 | 0 | 38.57 | 38.0 | 0.00 |
| emp-10k | B10k | dipper-r | 10000 | 0 | 23.97 | 23.0 | 0.02 |
| emp-10k | Bfull | sepp-prod | 10000 | 0 | 0.00 | 0.0 | 1.00 |
| emp-10k | Bfull | sepp-pinned | 10000 | 0 | 0.34 | 0.0 | 0.98 |
| emp-10k | Bfull | pplacer-scampp | 10000 | 0 | 8.69 | 0.0 | 0.68 |
| emp-10k | Bfull | epa-ng-scampp | 10000 | 0 | 9.30 | 0.0 | 0.55 |
| emp-10k | Bfull | bscampp | 10000 | 0 | 9.30 | 0.0 | 0.55 |
| emp-10k | Bfull | apples2 | 10000 | 0 | 18.19 | 12.0 | 0.16 |
| emp-10k | Bfull | dipper-m | 10000 | 0 | 78.08 | 76.0 | 0.00 |
| emp-10k | Bfull | dipper-m-gpu | 10000 | 0 | 69.89 | 65.0 | 0.00 |
<!-- /table -->

## 7. Recommendation

Basis: accuracy, runtime and memory only, read from the tables cited. Licenses are not a factor.

1. **Accuracy at the production reference size (Bfull\*).** No arm passes the pre-registered gate against `sepp-prod` (Table 5). The single pass (`pplacer-scampp`, V7-V9, n = 46) is an artefact of a tiny, noisy region (Section 3). `sepp-pinned`, which is SEPP with a newer pplacer, fails all nine regions, so the gate cannot be passed at Bfull\* by a near-identical tool at these sample sizes (Tables 3, 5), and this report does not call any arm equivalent to SEPP there. The effect sizes are the evidence: at full-length `bscampp` +2.26 [+1.44, +3.12], `epa-ng-scampp` +2.22 [+1.40, +3.08] and `pplacer-scampp` +1.95 [+1.23, +2.66] edges against a noise floor of -0.05 [-0.72, +0.60] (Table 5); in the other regions with at least 200 queries their intervals span 0 and are 2 to 6 edges wide, so a difference of about 2 edges can mostly be neither shown nor excluded (Table 5). `apples2` (+9.36), `dipper-m` (+27.65) and `dipper-m-gpu` (+55.20) are clearly worse at full-length (Table 5). No arm is picked as the best available replacement at Bfull\*.

2. **Accuracy at B10k\*.** `bscampp` and `epa-ng-scampp` pass in 5 of the 6 regions with at least 200 queries, `pplacer-scampp` in 4 of 6 (Table 3), and the simulated cross-check agrees for all three (Table 6). The V6-V8 fail of `bscampp` and `epa-ng-scampp` is a narrow miss (upper limits +0.42 and +0.38 against 0.25; Table 4) in a region where the `hmmalign` alignment costs the likelihood arms about 0.5 to 0.8 edges (Table 8). V1-V2, V1-V3 and V7-V9 (n = 153, 144, 46) do not decide anything. `apples2` fails all nine gg2 regions even with truth seeded by APPLES-2, so its accuracy counts for nothing here (Table 4). The DIPPER arms fail every gg2 region at both sizes and `dipper-r` fails at B10k\* (Tables 4, 5). The controls `epa-ng` (4 of 6) and `pplacer` (3 of 6) pass in no more regions than the SCAMPP-family arms (Table 3), and `sepp-pinned` is the version-drift control, not a challenger.

3. **Cost.** Among the challengers that clear the B10k\* gate, `bscampp` is the cheapest (the `epa-ng` control is cheaper still in some cells, Tables 13 to 20): on the same CPU it takes 0.15 to 0.23 of `sepp-prod`'s end-to-end time and 0.20 (full-length) to 0.80 of its memory at B10k\* (Tables 13, 15, 17, 19). `epa-ng-scampp` and `pplacer-scampp` are not reliably cheaper than SEPP in time or memory (Section 5). At Bfull\* the timing comparison with SEPP is confounded by hardware (SEPP on Xeon 6242, challengers on EPYC 7302; Table 2); `bscampp`'s end-to-end time is 0.15 to 0.73 of SEPP's and its memory 0.08 to 0.33 (Tables 14, 16, 18, 20), a range and not a hardware-normalised claim. Memory compares across machines: 22.8 GB against 277 GB at full-length (Table 14). `apples2` cannot be ranked end to end because its conversion fails.

4. **GPU.** `dipper-m-gpu` is among the cheapest placers (0.06 to 0.22 of SEPP end to end; Tables 13 to 20) but fails the accuracy gate, so it is reported beside the CPU ranking and does not compete in it.

5. **Conclusion.** This benchmark does not show that any arm is as accurate as production SEPP at Bfull\*, and at full-length the SCAMPP-family arms are about 2 edges worse on average there. It therefore does not support replacing production SEPP at the production reference size on accuracy. On B10k\* evidence and cost only (not as a Bfull\* pick), `bscampp` is the candidate to evaluate further, because it clears the gate at B10k\* in as many regions as any other challenger and is the cheapest of them. Whether about 2 edges of extra mean error at Bfull\* full-length (median 1.5 against 1.0; Table 7) is acceptable for the lab's use is a decision for the user that this report does not make.

## 8. Caveats

- **Model**: the shipped GTRCAT model has alpha = 1.0 and is used unchanged by every likelihood arm, as production SEPP does. It lowers absolute accuracy for all likelihood arms; the paired verdict against `sepp-prod` is unaffected (Tables 4, 5).
- **Pruned Bfull\***: it is the shipped reference minus the held-out queries (about 3% smaller), not the production reference itself. Bfull\* and B10k\* are never pooled (Tables 4, 5).
- **Truth bias**: GG2 truth was seeded by uDance (APPLES-2), so `apples2` gg2 results are optimistic; `apples2` still fails every gg2 region at B10k\* (Table 4). The non-uDance cross-check `sim10k` exists at B10k\* only (Table 6); there is no simulated check at Bfull\*, so Bfull\* results rest on gg2 truth alone.
- **DIPPER**: placement is joint (queries are placed together) and is scored through the Newick adapter; the adapter contracts zero-length backbone edges so that CPU `dipper-m` could be scored at Bfull\*. `dipper-m-gpu` is a different build (upstream main, not 0.1.5) and its accuracy differs from `dipper-m`, for example +55.20 against +27.65 at Bfull\* full-length (Table 5) and +2.45 against -0.38 oracle-minus-`hmmalign` at V3-V4 (Table 8). Both fail the gate everywhere.
- **hmmalign arm**: challengers use the deployable single-profile alignment while SEPP aligns itself; the alignment costs the likelihood arms up to about 1.3 edges in regions with at least 200 queries (Table 8, B10k\* only). The oracle alignment was not run at Bfull\*.
- **Short window**: V4-150 fragments fit several edges equally; SEPP and the likelihood arms all place only 0.17 of queries at distance 0 at B10k\* (Table 7), hence the paired comparison.
- **mem_source**: every placement row used in the cost tables has `mem_source=cgroup`; the 99 `psutil` rows are align and convert phases (Table 1), which Tables 13 to 20 do not use for peak memory.
- **Hardware**: SEPP at Bfull\* ran on Xeon 6242, challengers on EPYC 7302 (Table 2); B10k\* cells are all EPYC. Wall times against SEPP at Bfull\* are not hardware-normalised (Section 5).
- **DNFs are data**: 30 convert DNFs and no OOM or timeout (Table 1). 29 are `apples2` `guppy tog` errors (allow-listed); one is `pplacer` at gg2 Bfull\* V4-150, which placed in 2.6 h with 74 GB but whose likelihoods underflowed, so 0 of 879 queries are scored (`unplaced` in Table 5; `guppy tog` cannot read `-nan`; `docs/tools/pplacer.md`). That is the one capped `pplacer` Bfull\* cell.
- **Coverage**: `epa-ng`, `dipper-r` and `pplacer` are B10k\*-only; the oracle alignment and `sim10k` exist at B10k\* only (Tables 3, 6, 8).
- **Small regions**: V1-V2 (153 queries), V1-V3 (144) and V7-V9 (46) have intervals several edges wide and are never used to decide (Tables 4, 5).
- **Power at Bfull\***: the 0.25-edge margin is far below the width of the Bfull\* intervals (Table 5); the gate cannot show equivalence there.

## 9. Licenses

Flagged only, from `LICENSES.md`; no analysis, and licenses played no part in any ranking or in the recommendation. Copyleft engines to be reviewed later: GPL-3 (SEPP, pplacer, guppy, taxtastic, APPLES-2, TreeCluster) and AGPL-3 (EPA-ng, and therefore `epa-ng-scampp`, `bscampp` and `epa-ng`). `dipper-m` and `dipper-r` are MIT; the CUDA runtime in the GPU image is unverified (Table 22).

<!-- table:licenses -->
**Table 22.** Licenses (flag only; from LICENSES.md).

| Arm / tool | Version | Own license | Engine(s) executed | Copyleft flag |
|---|---|---|---|---|
| `sepp-prod` (q2-fragment-insertion) | 2026.7.0 | BSD-3 | SEPP 4.5.6 (GPL-3), pplacer 1.1.alpha19 (GPL-3), HMMER 3.4 (BSD-3), `guppy tog` (GPL-3), OpenJDK 25.0.2 | GPL-3 engines |
| `sepp-pinned` | SEPP 4.5.6 | GPL-3 | pplacer 1.1.alpha23 (GPL-3), HMMER 3.4 (BSD-3), OpenJDK 25.0.2 | GPL-3 |
| `pplacer` | 1.1.alpha23 | GPL-3 | itself | GPL-3 |
| `pplacer-scampp` | bscampp 1.0.8 | MIT (flag classifier says GPL) | pplacer (GPL-3), taxit from taxtastic (GPL-3), bundled `tools/hamming_distance/*` (MIT) | GPL-3 engines |
| `epa-ng-scampp` | bscampp 1.0.8 | MIT | EPA-ng 0.3.8 (AGPL-3), bundled `tools/hamming_distance/*` (MIT) | AGPL-3 engine |
| `bscampp` | bscampp 1.0.8 | MIT | EPA-ng 0.3.8 (AGPL-3), bundled `tools/hamming_distance/*` (MIT) | AGPL-3 engine |
| `epa-ng` (control) | 0.3.8 | AGPL-3 | itself | AGPL-3 |
| `apples2` | 2.0.11 | GPL-3 | itself, TreeCluster 1.0.5 (GPL-3); bundled FastTree is not run when `-D` is passed; ME branch lengths from pinned FastTree 2.2.0 (`FastTreeMP`) in the `ref-prep` phase | GPL-3 |
| `dipper-m` | 0.1.5 | MIT | itself | none flagged |
| `dipper-m-gpu` | image `swalia14/dipper@sha256:0b1fba6eef8651a7c9caaed5d07dbc0b57703824f40d6d0e1ad930fbc8d08fa2` (built from upstream main) | MIT | DIPPER itself; CUDA runtime in the image [UNVERIFIED] | none flagged (CUDA runtime license [UNVERIFIED]) |
| `dipper-r` | 0.1.5 | MIT | itself | none flagged |
| HMMER | 3.4 | BSD-3 | itself | none flagged |
| FastTree | 2.2.0 | source header says GPL-2+, repo `LICENSE` says GPL-3, bioconda says GPL-2.0-or-later | itself (`FastTree`, `FastTreeMP`) | GPL |
| guppy (from pplacer) | 1.1.alpha23 | GPL-3 | itself | GPL-3 |
| taxtastic (`taxit`) | 1.0.1 | GPL-3 | itself | GPL-3 |
| TreeCluster | 1.0.5 | GPL-3 | itself | GPL-3 |
| raxml-ng | 2.0.3 | AGPL-3 | itself (CI model-parse check only; not in any arm's run path) | AGPL-3 |
<!-- /table -->
