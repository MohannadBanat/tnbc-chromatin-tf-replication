# H3 (CRISPR dependency) results, 2026-11-06

Pre-specified in plan amendment 01 (SHA-256 60ace3eb...f461ff). Groups frozen before outcomes (SHA-256 75270bf4...cbb94).
Data: DepMap 23Q2 (Chronos), Project Score Chronos, DEMETER2; all files MD5-verified.

| Analysis | Result |
| --- | --- |
| H3 primary: HC-TF enrichment among TNBC dependencies (25 lines) | 5/93 (5.4%) vs 7.0% matched null (95% 2.2-11.8%); fold 0.77; p = 0.52 |
| H3a excluding common essentials | 0/87; all 5 hits (MED30, MYC, ILF2, HMGA1, SRF) common-essential; no TNBC-selective dependency at FDR <= 0.05 |
| H3a nominal TNBC < luminal | SOX9, E2F3 (p = 0.009 each, FDR 0.28) |
| H3b activity vs dependency within TNBC | median rho -0.023 vs size-matched null -0.020 (95% -0.072 to 0.039); p = 0.93; sign test p = 0.42 |
| H3b activity vs mRNA | mRNA tracks dependency slightly better: paired diff +0.064, Wilcoxon p = 0.027 |
| H3c Project Score (14 TNBC lines) | 7/93 vs 7.2%; p = 1.0; H3a 0/82; H3b median rho 0.00, p = 1.0 |
| H3c screen agreement (14 shared lines) | per-TF median effect rho = 0.82 (HC-TFs), 0.83 (all TFs) |
| H3d RNAi (DEMETER2, 32 TNBC lines) | 2/88 (MYBL2, PLAGL1); CRISPR-RNAi rho = 0.30 |
| H3e confounders | ZIC4 (only FDR-significant H3b hit): not a dependency (median +0.06), technical R2 = 0.08, not replicated in Project Score (rho -0.21, p = 0.47) |
| H3f subtypes | non-common-essential dependencies: FOXM1 (basal A, 7 lines); FOSL1, WWTR1 (basal unspecified, 3 lines) |

**Decision rule (amendment 01):** H3 is not positive, so working title reverts to #1.
