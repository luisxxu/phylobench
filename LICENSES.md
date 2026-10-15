# Licenses

Licenses are flagged only. This is not legal analysis and is not used for any scope or recommendation decision (design.md, principle 1). Wrappers never edit this file.

"Engine(s) executed" is what actually runs, not just what is installed.

| Arm / tool | Version | Own license | Engine(s) executed | Copyleft flag | Source at pinned tag |
|---|---|---|---|---|---|
| `sepp-prod` (q2-fragment-insertion) | 2026.7.0 | BSD-3 | SEPP 4.5.6 (GPL-3), pplacer 1.1.alpha19 (GPL-3), HMMER 3.4 (BSD-3), `guppy tog` (GPL-3), OpenJDK 25.0.2 | GPL-3 engines | https://github.com/qiime2/q2-fragment-insertion/tree/2026.7.0 |
| `sepp-pinned` | SEPP 4.5.6 | GPL-3 | pplacer 1.1.alpha23 (GPL-3), HMMER 3.4 (BSD-3), OpenJDK 25.0.2 | GPL-3 | https://github.com/smirarab/sepp/tree/v4.5.6 |
| `pplacer` | 1.1.alpha23 | GPL-3 | itself | GPL-3 | https://github.com/matsen/pplacer/tree/v1.1.alpha23 |
| `pplacer-scampp` | bscampp 1.0.8 | MIT (flag classifier says GPL) | pplacer (GPL-3), taxit from taxtastic (GPL-3), bundled `tools/hamming_distance/*` (MIT) | GPL-3 engines | https://github.com/ewedell/BSCAMPP/tree/v1.0.8 |
| `epa-ng-scampp` | bscampp 1.0.8 | MIT | EPA-ng 0.3.8 (AGPL-3), taxit (GPL-3), bundled `tools/hamming_distance/*` (MIT) | AGPL-3 engine | https://github.com/ewedell/BSCAMPP/tree/v1.0.8 |
| `bscampp` | bscampp 1.0.8 | MIT | EPA-ng 0.3.8 (AGPL-3), taxit (GPL-3), bundled `tools/hamming_distance/*` (MIT) | AGPL-3 engine | https://github.com/ewedell/BSCAMPP/tree/v1.0.8 |
| `epa-ng` (control) | 0.3.8 | AGPL-3 | itself | AGPL-3 | https://github.com/pierrebarbera/epa-ng/tree/v0.3.8 |
| `apples2` | 2.0.11 | GPL-3 | itself, TreeCluster 1.0.5 (GPL-3); bundled FastTree is not run when `-D` is passed | GPL-3 | https://github.com/balabanmetin/apples/tree/v2.0.11 |
| `dipper-m` | 0.1.5 | MIT | itself | none flagged | https://github.com/TurakhiaLab/DIPPER/tree/v0.1.5 |
| `dipper-r` | 0.1.5 | MIT | itself | none flagged | https://github.com/TurakhiaLab/DIPPER/tree/v0.1.5 |
| HMMER | 3.4 | BSD-3 | itself | none flagged | https://github.com/EddyRivasLab/hmmer/tree/hmmer-3.4 |
| FastTree | 2.2.0 | source header says GPL-2+, repo `LICENSE` says GPL-3, bioconda says GPL-2.0-or-later | itself | GPL | https://github.com/morgannprice/fasttree/tree/v2.2.0 |
| guppy (from pplacer) | 1.1.alpha23 | GPL-3 | itself | GPL-3 | https://github.com/matsen/pplacer/tree/v1.1.alpha23 |
| taxtastic (`taxit`) | 1.0.1 | GPL-3 | itself | GPL-3 | https://github.com/fhcrc/taxtastic/tree/v1.0.1 |
| TreeCluster | 1.0.5 | GPL-3 | itself | GPL-3 | https://github.com/niemasd/TreeCluster/tree/1.0.5 |
| raxml-ng (T6 pins and adds its env) | pinned in T6 | AGPL-3 | itself | AGPL-3 | https://github.com/amkozlov/raxml-ng |

## Bundled-binary exception

BSCAMPP/SCAMPP use the pinned external `pplacer`, `epa-ng`, and `taxit` (never their bundled copies), but the bundled `bscampp/tools/hamming_distance/{fragment_hamming,hamming,homology}` executables do run.
