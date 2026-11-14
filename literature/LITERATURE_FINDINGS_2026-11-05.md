# Literature search: findings (run 2026-11-05)

**Scope:** 11 PubMed searches (S1–S10), 989 hits, 984 unique records, title-screened (abstracts for borderline titles). Plus citation chasing (OpenAlex, Semantic Scholar) and 12 hand-picked methods papers. **100 records kept** (49 A, 44 B, 7 C). All 98 kept DOIs resolve (97 Crossref, 1 DataCite). Full log: `literature_log.xlsx`.

Screening decisions are a first pass by Claude and need a human check before the reference list is final.

## Novelty checks

| Question | Finding | Status |
|---|---|---|
| Has anyone replicated or critiqued Bahl et al. 2026? | 0 citations in OpenAlex and Semantic Scholar | Checked; **Google Scholar not checked** (no API). Re-check before submission |
| Has a TNBC TF set been tested systematically against CRISPR dependencies (H3)? | S2 found only single-gene studies (ID4, NFIB, SNAI1, FOXC1, MYC, ZBTB7A). Closest broader work: Marcotte 2016 (shRNA, breast lines), Lehmann 2021 (TNBC subtype vulnerabilities) | No systematic TF-set test found in PubMed. Not proof of absence: bioRxiv/medRxiv not searched yet |
| Has the Lund group published on GSE337499/500? | None of the 31 citers of Aine 2025 analyses these cell lines (by title) | Checked by title only |
| Cross-lab reproducibility of TNBC cell-line ATAC (H1)? | No direct study found. Harris 2023 (TNBC lines, by donor ancestry) is a candidate third ATAC dataset; Cell Rep Methods 2026 (breast cell-line reference samples) is the closest benchmark | Line overlap with Harris 2023 still to check |

## Records that change the plan

1. **Phipson & Smyth 2010, "Permutation P-values should never be zero."** Bahl et al. report p = 0 for 58 drug–TF pairs. This supports the PDX/C-index critique directly.
2. **Schröder et al. 2011 (survcomp) and Smirnov et al. 2016 (PharmacoGx).** The likely source of Bahl's C-index and its confidence intervals. Their Code Ocean capsule (verified, DataCite) should confirm which estimator was used.
3. **Mer et al. 2019, Cancer Res (Xeva).** The likely source of the "BestAvgResponse" PDX metric. Cite it in the PDX re-test methods.
4. **Cancer Res Commun 2026: multiple-comparison tests in small-sample PDX studies; NCI PDXNet consensus (2024); Neuro-Oncology 2021 on PDX design.** Current standards against which an n = 11 claim can be judged.
5. **Harris et al. 2023 (Cancer Res Commun, doi:10.1158/2767-9764.CRC-23-0236).** Possible third lab for H1. If it shares lines with Bahl and Lund, the cross-lab test gets much stronger.
6. **Pharmacogenomic consistency literature (Haibe-Kains 2013; the Nature 2016 replies; Ben-David 2018; RP:CB replication of GDSC markers).** Frames H4: drug–biomarker links often fail across screens, so independent replication is the standard.
7. **Competing mTOR biomarker: TOMM40/FADS2 ratio (Cancer Med 2026).** Must be discussed alongside NFE2L3/CBFB.

## Gaps still open

- bioRxiv/medRxiv search (last 24 months): not run yet.
- Google Scholar "Cited by" for Bahl et al.: not checkable by API; run by hand.
- S8a (C-index, small samples) returned almost nothing useful (1 of 98). Methods papers were added by hand instead; a Google Scholar search on "concordance index" + "small sample" could fill this.
- Recent Biomolecules/IJMS TNBC papers (journal-specific citations): only 2 found so far (Cells 2024, IJMS 2024).
