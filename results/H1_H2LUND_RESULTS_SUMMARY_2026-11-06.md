# H1 and H2 (Lund) results, 2026-11-06 (scripts 11b, 14; verified against results/ files)

**H2, Lund RNA-seq (GSE337499), 6 shared lines:** median per-TF rho 0.83 (bootstrap 95% CI 0.60-1.00); exact permutation p = 0.0014 (719 permutations); 73/93 HC-TFs rho >= 0.5. CBFB 0.94, NFE2L3 0.43. Exploratory line identity: released profile best matched by own Lund line for 3/6.

**H1, cross-lab chromatin** (Bahl GSE323120 peaks, 48 samples / 24 lines; Lund GSE337500 RPKM bigWigs, 16 samples / 8 lines; hg38). Shared lines: BT549, HCC1937, MDA-MB-231, MDA-MB-468.
- Consensus peaks (504,081; called in >= 2 Bahl samples): pre-specified criterion met 1/4 (MDA-MB-468 only, rho 0.69). Same-line cross-lab rho 0.39-0.69; within-lab different-line maxima often higher (lab effect).
- HC-TF promoters (TSS +/- 1 kb): criterion met 3/4 (BT549 0.84 vs 0.83; HCC1937 0.80 vs 0.78; MDA-MB-468 0.83 vs 0.77; MDA-MB-231 fails 0.72 vs 0.75). Margins 0.01-0.06.
- Exploratory identity: same line ranked first among 24 Bahl lines for 3/4 (raw) and 4/4 (centred within lab; rho 0.33-0.48 vs <= 0.29 for any other).
- Replicate agreement: Lund 0.74-0.88; Bahl 0.68-0.89.
- Verdict: H1 partly supported. Caveat: different deposited data types (peaks vs coverage).
