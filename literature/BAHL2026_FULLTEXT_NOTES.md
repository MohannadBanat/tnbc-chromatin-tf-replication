# Bahl et al. 2026 (Nat Commun 17:9897): full-text review notes

Source: published PDF supplied by the user on 2026-11-06 (16 pages). Page numbers refer to the journal pagination.

## Points that bear on our study

1. **The PDX C-index method is not stated.** Methods ("In vivo drug sensitivity", p. 12) say only that "the CI was also calculated between the response to everolimus and the activity scores of TFs in these PDXs", with BH correction. No estimator, no confidence-interval method and no software are named. For the cell lines they name PharmacoGx (ref. 38), but not for the PDXs. This explains why our recomputed C-indices differ slightly from theirs, and it makes their reported p-values unverifiable from the paper alone. → Inspect the Code Ocean capsule.
2. **No sample-size justification.** "No statistical method was used to predetermine sample size" (p. 12). The everolimus claim rests on n = 11 PDXs, from a different PDX cohort (ref. 64, Savage et al. 2020), not their own.
3. **The "three independent" drug datasets share most cell lines.** GRAY (24 lines), GDSC2 (19) and CTRPv2 (15) overlap: 13 lines are in all three, 14 of the 15 CTRPv2 lines are in GDSC2, and 17 of the 19 GDSC2 lines are in GRAY (Supp. Data 7). "Concordant in ≥2 datasets" therefore largely means the same cell lines re-screened in different assays, not independent biological samples. The Discussion calls this "robust biology rather than dataset-specific artifacts" (p. 7).
   - **Implication for our H4:** PRISM screens many of the same TNBC lines. Our replication is assay-independent, not sample-independent, and the manuscript must say so. Add to `DEVIATIONS.md` as an interpretation note (the analysis itself is unchanged).
4. **60 vs 94 TFs.** Methods say the drug CI used "activity scores of 60 HC-TFs" (p. 12). Results and figures describe 94 HC-TFs (Fig. 5A–C legend: "activity scores of 94 HC-TFs"). Supp. Data 7 contains 76 distinct TFs among its significant pairs. The number of TFs tested (which sets the multiple-testing burden) is unclear.
5. **Everolimus cell-line support comes from one dataset.** "Everolimus was also tested in TNBC cell lines in only the GRAY dataset" (p. 7), with NFE2L3 concordant there (Supp. Fig. 7D). So the in-vitro backing for the headline mTOR claim is a single screen.
6. **An unclear Methods sentence** (p. 12): PDX activity scores were "based on RNA-seq data recently released by Savage et al. 2020 as a reference set of position weight matrices", which mixes up the RNA-seq source with the motif reference. This is a likely editing error; it should be checked in the capsule, not cited as evidence.
7. **Claims to quote accurately:** "these results provide proof-of-concept that HC-TF activity, particularly NFE2L3 and CBFB, can stratify mTOR inhibitor response in preclinical models" (p. 7).

## Not yet read in detail
Supplementary Information PDF (Supp. Figs 1–7), the Reporting Summary and the Transparent Peer Review file. The peer-review file may show whether reviewers raised the small-n issue.
