# Chromatin-informed TF biomarkers in TNBC: an independent bioinformatics replication

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210657.svg)](https://doi.org/10.5281/zenodo.24210657)

Release v1.0 (the version submitted with the manuscript) is archived at Zenodo: https://doi.org/10.5281/zenodo.24210657

Code, pre-registered analysis plan and results for:

> Arabiat A., Abu Kiwan A., Banat M., Abbas M., Abu Owida H. *Integrative Bioinformatics Analysis of Chromatin-Informed Transcription-Factor Biomarkers in Triple-Negative Breast Cancer: An Independent Replication Study.* Submitted to *Biomolecules* (MDPI), 2026.

The study re-tests the 94 high-confidence TNBC transcription factors (HC-TFs) of Bahl et al. (*Nat. Commun.* 2026, doi:10.1038/s41467-026-76385-8) using only public data.

## Analysis plan
- `ANALYSIS_PLAN_2026-11-05.md` (SHA-256 in `.sha256`), committed before any outcome data were downloaded.
- `ANALYSIS_PLAN_AMENDMENT_01_2026-11-06.md`: CRISPR-dependency hypotheses (H3a–g) and the title decision rule.
- `groups/`: cell-line groups frozen from metadata before outcome data were read.
- `DEVIATIONS.md`: every departure from the plan, dated.

## Hypotheses → scripts → outputs (`results/`)
| | Question | Script(s) |
|---|---|---|
| PDX re-test (exploratory) | Do the everolimus associations survive exact permutation testing? | `01` |
| H1 | Is chromatin accessibility of shared lines reproducible across laboratories? | `14`, `14b` (sensitivity) |
| H2 | Is HC-TF activity reproducible in independent RNA-seq? | `11` (DepMap), `11b` (Lund) |
| H3 | Are HC-TFs CRISPR/RNAi dependencies of TNBC cells? | `04` (groups), `05`, `06`, `07` (Project Score), `08`, `09` |
| H4 | Do the reported drug–TF pairs replicate in PRISM and GDSC1? | `10`, `10b` (PRISM; null, 10,000 permutations), `16` (GDSC1; run `16 <root> 1` and `16 <root> 2` for the two null halves, then `16 <root> combine`) |
| Figures / tables | | `12`, `12b`, `13` |

Literature search: `LITERATURE_SEARCH_PLAN.md`, `literature/`, scripts `02`, `03`.

## Data (not redistributed; downloaded by DOI/accession)
DepMap 23Q2 (doi:10.6084/m9.figshare.22765112.v4), Project Score Chronos (doi:10.6084/m9.figshare.14461980.v2), DEMETER2 (doi:10.6084/m9.figshare.9170975.v1), PRISM Repurposing 20Q2 (doi:10.6084/m9.figshare.20564034.v1), GDSC1 fitted dose–response data (doi:10.6084/m9.figshare.19345397; file and checksum in `data/gdsc_raw/`), GEO GSE323120, GSE337499, GSE337500; Bahl et al. Supplementary Data (in `data/bahl2026/`, CC BY 4.0). Expected file names and checksums: `data/depmap_raw/manifest*.tsv` and `SHA256SUMS`.

## Reproducing
Python 3.13 with NumPy, pandas, SciPy, statsmodels, decoupler 2.3.1, pyBigWig and matplotlib. Each script takes the project root as its first argument, e.g. `python scripts/05_h3_primary_h3a.py .`. Seeds are fixed (20261105, 20261106, 202611101).

## Licence
Code: MIT. Bahl et al. supplementary files remain under their original CC BY 4.0 licence.
