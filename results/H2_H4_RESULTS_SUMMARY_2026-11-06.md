# H2 and H4 results (2026-11-06), verified against results/ files

**H2 (activity reproducibility, DepMap 23Q2 RNA-seq, 49 shared lines):** median per-TF Spearman rho 0.80 (bootstrap 95% CI 0.77-0.84) across 93 HC-TFs; line-permutation null median -0.003, p = 0.001. 86/93 rho >= 0.5, none < 0.2. Other TNBC-specific TFs median 0.74. NFE2L3 0.88, CBFB 0.87. Lund RNA-seq (GSE337499) part pending: GEO blocked.

**H4 (PRISM Repurposing 20Q2 secondary, doi:10.6084/m9.figshare.20564034.v1):** 176 pairs reported in >= 2 datasets, 148 testable, up to 16 TNBC lines. 3 replicated (2.0%; 95% CI 0.4-5.8%): afatinib-GRHL3, dasatinib-HMGB3, tozasertib-LMO1. Random-pair null 1.9% (0-4.1%), p = 0.55.
mTOR pre-specified: 28 compounds x NFE2L3/CBFB = 56 tests, min FDR 0.375. Everolimus: NFE2L3 C 0.475 (p 0.83), CBFB C 0.533 (p 0.75), n 16. Sirolimus-CBFB C 0.80, p 0.0067, n 12 (FDR 0.375); other rapalogs C 0.53-0.64, p >= 0.28.

**H1:** not run; GEO (ftp.ncbi.nlm.nih.gov) not on the workspace allowlist.
