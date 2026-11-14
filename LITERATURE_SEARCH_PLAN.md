# Literature search plan: TNBC TF-dependency replication study

Target journal: Biomolecules (MDPI). Created 2026-11-05.
Search window: 2010 – present (no limit for the methods blocks C and F).
Databases: PubMed, Scopus (or Web of Science), Google Scholar (first 200 results), bioRxiv/medRxiv (last 24 months).
For every search, record: database, date run, exact string, number of hits, number kept.

---

## 1. Concept blocks (combine with AND between blocks, OR within)

| Block | Concept | Keywords / synonyms |
|---|---|---|
| **A** | Disease | "triple-negative breast cancer" OR TNBC OR "basal-like breast cancer" OR "ER-negative breast cancer" |
| **B** | Transcription-factor activity | "transcription factor activity" OR "master regulator" OR regulon OR VIPER OR ARACNe OR "protein activity inference" OR "gene regulatory network" OR "transcriptional network" |
| **C** | Chromatin | "chromatin accessibility" OR ATAC-seq OR "open chromatin" OR "cis-regulatory element" OR enhancer OR chromVAR OR "motif enrichment" |
| **D** | Genetic dependency | "CRISPR screen" OR "gene dependency" OR "essential genes" OR DepMap OR "Cancer Dependency Map" OR Chronos OR CERES OR "synthetic lethality" |
| **E** | Drug response | pharmacogenomic* OR "drug sensitivity" OR GDSC OR CTRP OR PRISM OR "cell line drug screen" OR "drug response prediction" |
| **F** | Reproducibility | reproducib* OR replicat* OR "independent validation" OR concordance OR "inter-laboratory" OR "batch effect" OR "re-analysis" |
| **G** | Specific targets | everolimus OR "mTOR inhibitor" OR NFE2L3 OR CBFB OR rapalog |
| **H** | Small-n statistics | "concordance index" OR "C-index" OR "permutation test" OR "small sample" OR "patient-derived xenograft" |

---

## 2. Searches, one per claim in the paper

**S1 – Prior work: TF activity + chromatin in TNBC** (what we build on; novelty check)
PubMed: `("triple-negative breast cancer"[tiab] OR TNBC[tiab] OR "basal-like"[tiab]) AND ("transcription factor activity"[tiab] OR "master regulator"[tiab] OR regulon[tiab] OR VIPER[tiab]) AND ("chromatin accessibility"[tiab] OR ATAC-seq[tiab] OR enhancer[tiab])`

**S2 – Have TNBC TFs been tested as CRISPR dependencies?** (novelty of H3)
PubMed: `(TNBC[tiab] OR "triple-negative"[tiab] OR "basal-like"[tiab]) AND ("transcription factor"[tiab]) AND (CRISPR[tiab] OR DepMap[tiab] OR "gene dependency"[tiab] OR "essential gene*"[tiab])`

**S3 – TF activity as a drug-response biomarker** (context for H4)
PubMed: `("transcription factor activity"[tiab] OR "master regulator"[tiab] OR VIPER[tiab] OR "protein activity"[tiab]) AND ("drug sensitivity"[tiab] OR pharmacogenomic*[tiab] OR GDSC[tiab] OR CTRP[tiab] OR PRISM[tiab])`

**S4 – Pharmacogenomic dataset agreement** (why independent screens matter)
PubMed: `(GDSC[tiab] OR CTRP[tiab] OR CCLE[tiab] OR PRISM[tiab] OR "cancer cell line"[tiab]) AND (consistency[tiab] OR concordance[tiab] OR reproducib*[tiab] OR inconsisten*[tiab]) AND ("drug response"[tiab] OR "drug sensitivity"[tiab])`
Key known debate: Haibe-Kains et al. 2013 (Nature) and the responses that followed.

**S5 – Cross-laboratory ATAC-seq reproducibility** (H1)
PubMed: `(ATAC-seq[tiab] OR "chromatin accessibility"[tiab]) AND (reproducib*[tiab] OR "batch effect*"[tiab] OR "inter-laboratory"[tiab] OR "cross-laboratory"[tiab] OR benchmark*[tiab])`

**S6 – Robustness of regulatory-network (VIPER/ARACNe) inference** (H2)
PubMed: `(VIPER[tiab] OR ARACNe[tiab] OR DoRothEA[tiab] OR decoupleR[tiab] OR regulon*[tiab]) AND (benchmark*[tiab] OR reproducib*[tiab]) AND ("transcription factor*"[tiab])`

**S7 – mTOR inhibitors in TNBC: biomarkers and trials** (context for the everolimus re-test)
PubMed: `(everolimus[tiab] OR "mTOR inhibitor*"[tiab] OR rapalog*[tiab]) AND (TNBC[tiab] OR "triple-negative"[tiab]) AND (biomarker*[tiab] OR predict*[tiab] OR "clinical trial"[pt])`

**S8 – Statistics for small n / concordance index** (methods justification)
PubMed: `("concordance index"[tiab] OR "C-index"[tiab]) AND ("small sample*"[tiab] OR permutation[tiab])`
Also: `("patient-derived xenograft"[tiab]) AND ("sample size"[tiab] OR "statistical power"[tiab] OR "study design"[tiab])`

**S9 – Replication and research-waste context in cancer biology** (Introduction framing)
PubMed: `("Reproducibility Project: Cancer Biology"[tiab]) OR ((replication[tiab] OR reproducib*[tiab]) AND "cancer biology"[tiab] AND (preclinical[tiab]))`

**S10 – TNBC cell-line models: what they represent** (Discussion)
PubMed: `(TNBC[tiab] OR "triple-negative"[tiab]) AND ("cell line*"[ti]) AND (subtype*[tiab] OR Lehmann[tiab] OR "molecular classification"[tiab])`

---

## 2b. PubMed hit counts (all strings tested on 2026-11-05)

| Search | Hits | Note |
|---|---|---|
| S1 TF activity + chromatin, TNBC | 10 | Screen all. Contains Bahl 2026 and Harris 2023 (TNBC cell-line ATAC by donor ancestry), a possible third ATAC dataset for H1 |
| S2 TNBC TFs as CRISPR dependencies | 29 | Screen all. Only single-gene studies (ID4, FOXC1, SNAI1, ZBTB7A…); none tests a TF set systematically against DepMap |
| S3 TF activity → drug response | 63 | Screen all |
| S4 Pharmacogenomic consistency | 60 | Screen all |
| S5 ATAC-seq reproducibility | 364 | Screen titles; keep benchmarks and cross-lab studies |
| S6 Regulon-method benchmarks | 38 | Narrowed from 823 |
| S7 mTOR inhibitors in TNBC | 56 | Screen all |
| S8a C-index, small samples/permutation | 98 | Narrowed from 3,371 |
| S8b PDX sample size/power | 35 | Screen all |
| S9 Cancer-biology reproducibility | 92 | Keep the Reproducibility Project papers plus key commentary |
| S10 TNBC cell-line subtypes | 144 | Narrowed from 1,642 |

About 990 records before de-duplication; roughly 600 titles to screen once S5 is title-filtered.

## 3. Citation chasing (do these by hand; often more useful than keyword searches)

1. **Forward citations of Bahl et al. 2026** (doi:10.1038/s41467-026-76385-8): has anyone already replicated or criticised it? Use Google Scholar "Cited by" and Scopus. **Most important search; repeat it right before submission.**
2. Backward references of Bahl et al.: their refs 24–25 (ARACNe3/VIPER), 35–39 (GRAY, GDSC2, CTRPv2, PharmacoGx, ORCESTRA), 41 (PDX everolimus source).
3. Forward citations of Aine et al. 2025, Nat Commun (doi:10.1038/s41467-025-58158-x) and of the GSE337499/GSE337500 data (Lund lines). Checks whether the Lund group has published on these lines yet.
4. The DepMap Chronos paper (Dempster et al. 2021, Genome Biology) and the PRISM Repurposing paper (Corsello et al. 2020, Nature Cancer): methods citations.
5. Recent TNBC papers in Biomolecules and IJMS (search within the journal sites for "triple-negative" + "transcription factor"). Shows editors' expectations and gives citable journal-specific work.

---

## 4. Screening rules

- **Keep if** it (a) tests TF activity or TF dependency in breast cancer models, (b) assesses reproducibility of omics or pharmacogenomic findings, (c) is a primary method we use (VIPER, ARACNe, chromVAR, Chronos, PRISM, C-index), or (d) is a direct competitor (same question, any cancer type).
- **Drop if** it is a single-gene wet-lab study with no network/regulatory angle, a prognostic gene-signature paper without validation, or a review older than 5 years (unless it is the standard reference).
- Log decisions in `literature_log.xlsx`: ID | citation | DOI | search (S1–S10) | keep/drop | reason | which paper section.
- Verify every DOI resolves before it enters the manuscript.

## 5. What would change the plan

- A paper already testing Bahl's HC-TFs against DepMap, or re-testing their drug associations → reframe around whatever remains new.
- A Lund publication analysing GSE337499/500 with DepMap → cite it, and shift emphasis to the cross-lab comparison.
