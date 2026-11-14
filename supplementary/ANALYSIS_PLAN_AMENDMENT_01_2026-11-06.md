# Analysis plan amendment 01 — committed on 2026-11-06

Amends: internal analysis plan committed on 2026-11-05 (SHA-256 b69fd156912c863fc546c1550903ab719ce8ceee4660cad35e6cc92252e6cf37).

**Reason:** the authors chose CRISPR dependency (H3) as the paper's main axis (working title: "Independent Validation of Transcription-Factor Dependencies in Triple-Negative Breast Cancer: A Computational Replication Study of a Chromatin-Informed Regulatory Framework"). This amendment deepens H3. H1, H2 and H4 are unchanged and become supporting analyses.

**Blinding statement:** written before any DepMap, Project Score or RNAi file (gene effects, expression or model metadata) was downloaded or opened. The only DepMap-related information seen so far is from Bahl et al.'s supplementary files.

---

## H3 primary analysis (unchanged from the original plan)

The proportion of the 94 HC-TFs with median Chronos gene effect ≤ −0.5 across TNBC/basal-like lines, compared with 10,000 expression-decile-matched draws from the 1,996 TFs of the TCGA network. Two-sided empirical p.

## New pre-specified H3 analyses

**H3a — Separate selective from common-essential dependencies (key confound).**
Some HC-TFs may be essential in nearly every cell line (e.g., cell-cycle or general transcription regulators). Enrichment driven by them says nothing about TNBC.
- Flag HC-TFs on DepMap's common-essential list for the release used.
- Repeat the primary test excluding common essentials, with the null also excluding them.
- **TNBC selectivity:** gene effect in TNBC/basal-like lines vs (i) luminal breast lines and (ii) all other lineages. Wilcoxon test plus the difference in medians with a bootstrap 95% CI. BH-FDR ≤ 0.05 across HC-TFs.
- An HC-TF is a "TNBC-selective dependency" only if (a) its median effect in TNBC ≤ −0.5, (b) it is not common-essential, and (c) it passes the selectivity test against luminal lines.

**H3b — Does activity predict dependency? (the core test of the framework)**
Within TNBC/basal-like lines, per HC-TF: Spearman correlation between VIPER activity (Bahl's released regulons, applied to DepMap expression) and Chronos gene effect. A negative correlation means more active lines depend on the TF more.
- Report the distribution over HC-TFs and a one-sample sign test of median correlation < 0, against 1,000 size-matched regulon sets.
- Secondary: the same with the TF's own mRNA instead of activity. Activity should beat expression if the framework adds value; we report the paired difference.

**H3c — Replication in an independent CRISPR screen.**
Repeat the primary analysis, H3a and H3b in the Sanger Project Score CRISPR data, using the version processed with Chronos where DepMap provides it, otherwise the original scores. Report agreement of per-TF median gene effect between screens (Spearman, lines screened in both).

**H3d — Orthogonal perturbation modality.**
Repeat H3a's dependency calls with RNAi (DEMETER2). Report how many TNBC-selective CRISPR dependencies are also RNAi dependencies. Descriptive only; RNAi off-target effects are acknowledged.

**H3e — Confounders.**
- **Copy number:** for each selective dependency, test whether the gene effect is explained by the TF's copy number (CRISPR cutting toxicity in amplified regions). Partial Spearman adjusting for copy number.
- **Technical-only negative-control model** (from the original plan): growth rate, screen quality, library size. Flag associations this model explains as well as activity does.
- **Paralog buffering (descriptive):** for HC-TFs in families (e.g., TEAD2/TEAD4, E2F2/E2F3/E2F8), report whether paralog expression relates to dependency strength.

**H3f — Subtype structure.**
Assign lines to basal-like A, basal-like B (mesenchymal) or luminal-androgen-receptor using published cell-line annotations, defined before the gene-effect data are opened. Report per-subtype median effects for TNBC-selective dependencies. Descriptive; small groups get permutation p-values only.

**H3g — The two laboratories' cell lines.**
For the cell lines in GSE323120 and GSE337500 that DepMap screened, show HC-TF activity (from each laboratory's data) next to the dependency, as an illustrative figure. Not a hypothesis test.

## Decision rules for the title

- **Positive H3** (title #2 stands): the primary test shows enrichment at p < 0.05 **and** it survives excluding common essentials (H3a) **and** it replicates in direction in Project Score (H3c).
- Otherwise use title #1. Any partial result (e.g., enrichment entirely due to common essentials) is reported as the finding, not hidden.

## Unchanged

Unit of analysis = cell line; permutation p-values whenever fewer than 8 lines are in a test; BH-FDR within each analysis family; effect sizes with 95% bootstrap CIs; no pooling across screens; no predictive model fitted (nested CV not needed).
