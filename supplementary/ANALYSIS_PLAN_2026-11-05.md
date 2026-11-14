# Internal analysis plan committed on 2026-11-05

**Working title:** Are chromatin-informed transcription-factor dependencies in triple-negative breast cancer reproducible? An independent cross-laboratory, CRISPR-dependency and pharmacogenomic validation study

**Study type:** Retrospective computational re-analysis and independent validation study using public data only (no new samples).

**Target journal:** Biomolecules (MDPI), Bioinformatics and Systems Biology area; alternative: IJMS.

**Claims under test:** Bahl et al., *Nat Commun* 2026 (doi:10.1038/s41467-026-76385-8): (a) 94 high-confidence TNBC transcription factors (HC-TFs) are defined by VIPER activity plus chromatin accessibility; (b) HC-TF activity is associated with drug sensitivity across GRAY, GDSC2 and CTRPv2 (176 drug–TF pairs significant in ≥2 datasets); (c) NFE2L3/CBFB activity predicts everolimus sensitivity in 11 PDXs.

This plan is not registered on OSF. It is an internal plan, hashed (SHA-256) and stored in the project repository before the outcome data below are downloaded or examined.

---

## 0. What was examined before this plan (disclosed)

1. Bahl et al. Supplementary Data 1–10 (structure, HC-TF list, released VIPER network, released activity scores).
2. **Exploratory re-test of their everolimus PDX result** (script `01_recompute_everolimus_pdx.py`; results `01_everolimus_pdx_recomputed.csv`). Using their released activities and responses (n = 11 PDXs), a permutation test (200,000 permutations, seed 20261105) found 0/94 TFs with BH-FDR ≤ 0.05 (minimum unadjusted p = 0.041), versus 11/94 that they reported. This is reported as an exploratory reproducibility finding, not a pre-specified test.
3. Attempted recomputation of their 403 cell-line drug–TF concordance values from released data: cell-line counts matched, but the C-indices did not (median absolute difference 0.17). That indicates their drug analyses used activity scores not included in the supplements. **No conclusion is drawn from this.** Their CodeOcean capsule (doi:10.24433/CO.3288699.v1) will be inspected before any statement about those values.

Not yet downloaded or examined: GSE337499, GSE337500, GSE323120 peak/signal files; DepMap CRISPR (Chronos) gene effects; DepMap/CCLE expression; PRISM Repurposing; GDSC1.

## 1. Data

| Source | Use |
|---|---|
| GSE323120 (Bahl; ATAC-seq, 26 TNBC + 2 ER+ lines, 2 replicates) | Lab A chromatin |
| GSE337500 (Lund; ATAC-seq, 8 TNBC lines, 2 replicates) | Lab B chromatin |
| GSE337499 (Lund; RNA-seq, same 8 lines) | Independent expression for VIPER |
| Bahl Supplementary Data 5 (94-HC-TF VIPER network, 6,844 edges; cell-line activities) | Frozen regulons and reference activities |
| DepMap (latest public release at download; release ID and checksums recorded) | CRISPR gene effect (Chronos), expression |
| PRISM Repurposing (secondary screen) | Independent drug-response assay |
| GDSC1 | Secondary drug dataset, not used by Bahl et al. |
| GSE176078 (single-cell TNBC atlas) | Secondary: confirm HC-TF activity is tumour-cell-intrinsic |

Shared lines: BT-549, HCC1937, MDA-MB-468 (both labs, in Bahl's included peak set); MDA-MB-231 (in GSE323120 but not in Bahl's included peak set; analysed as a sensitivity case).

**Unit of analysis: the cell line.** Replicates are averaged before any between-line statistic.

## 2. Hypotheses and primary analyses

**H1 – Cross-laboratory chromatin reproducibility.** For each shared line, accessibility at HC-TF promoters (±1 kb TSS) and genome-wide consensus peaks is more similar between labs for the same line than between different lines within one lab.
- Processing: identical pipeline for both datasets from FASTQ/BAM where available (same aligner, same peak caller, IDR across replicates); otherwise from the deposited bigWig signal at a common peak set.
- Metric: Spearman correlation of log-CPM at consensus peaks; Jaccard of IDR peaks.
- Success criterion (fixed now): same-line cross-lab correlation is greater than every different-line within-lab correlation for ≥3 of 4 shared lines.
- Secondary: chromVAR motif deviation for HC-TFs with motifs (56/94); cross-lab Spearman per TF across shared lines.

**H2 – Reproducibility of HC-TF activity.** Apply Bahl's released 94-TF regulons unchanged (VIPER, same parameters as their methods) to (i) Lund RNA-seq (GSE337499) and (ii) DepMap expression. Compare with their released cell-line activities.
- Metric: per-TF Spearman across shared lines; median over TFs, with 95% bootstrap CI (resampling lines).
- Null: the same statistic for 94 regulons of matched size drawn from the full TCGA VIPER network (1,000 draws).

**H3 – Are HC-TFs genetic dependencies? (primary new question)** In DepMap TNBC/basal-like lines:
- Primary: the proportion of the 94 HC-TFs with median Chronos gene effect ≤ −0.5, compared with the matched null (TFs from the 1,996-TF network, matched on expression decile; 10,000 draws). Two-sided empirical p.
- Secondary: per-TF Spearman between VIPER activity and gene effect across lines (a negative value means a dependency where active), BH-FDR.
- Selectivity: the same tests in luminal lines; an HC-TF counts as "TNBC-selective" only if the effect differs (Wilcoxon, BH-FDR ≤ 0.05).

**H4 – Independent replication of drug–TF associations.** For the 176 drug–TF pairs significant in ≥2 of their datasets, re-test those whose drug is present in PRISM and/or GDSC1:
- Statistic: concordance index between activity and AUC/AAC (same direction convention as Bahl), with a permutation p-value (10,000 permutations; the cell line as unit).
- Replication = same sign and permutation p < 0.05.
- Primary outcome: replication rate with a Clopper–Pearson 95% CI, compared with the replication rate of 1,000 random drug–TF pair sets of equal size.
- Pre-specified: everolimus/rapalogs/other mTOR inhibitors × NFE2L3 and CBFB.

## 3. Controls and robustness

- **Technical-only negative control model:** for H3/H4, a model using only technical covariates (library size, peak count, FRiP/TSS enrichment, DepMap screen quality, growth rate) is fitted alongside. Any association that the technical-only model explains equally well is flagged.
- **Nested cross-validation:** used only if a predictive model is fitted (e.g., an HC-TF dependency classifier); outer 5-fold × 20 repeats, inner tuning. No predictive claim is made without it.
- **Platform split:** H2 and H4 are reported separately per expression source (Lund RNA-seq vs DepMap) and per drug screen; no pooling across screens.
- **Multiple testing:** BH-FDR within each hypothesis family.
- **Effect sizes:** all with 95% CIs (bootstrap over cell lines, 2,000 resamples) and bootstrap optimism correction for any fitted model.
- **Small-n rule:** with fewer than 8 lines in a test, only permutation p-values are reported; no asymptotic standard errors.

## 4. Outputs

- The pre-specified tables and figures for H1–H4, plus the exploratory PDX re-test (Section 0).
- Before submission, the authors of Bahl et al. are contacted with the PDX re-test so that any explanation (e.g., a different C-index definition) can be addressed in the manuscript.
- Code: GitHub, MohannadBanat organisation, MIT licence; repository private until submission. All accessions, DOIs and links checked before submission.
- Reporting follows the relevant parts of REMARK (marker evaluation) for H3/H4.

## 5. Deviations

Any deviation from this plan will be listed in a dated `DEVIATIONS.md` and reported in the manuscript.
