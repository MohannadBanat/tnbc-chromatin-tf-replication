"""Build the literature screening log from the PubMed records + citation chasing.
Screening level: titles for all 984 records; abstracts consulted for borderline titles.
Decisions below are the screener's (Claude, 2026-11-05) and must be checked by a human."""
import sys, pandas as pd
from pathlib import Path

lit = Path(sys.argv[1])
df = pd.read_csv(lit / 'pubmed_records_raw_2026-11-05.csv', dtype=str).fillna('')

# pmid: (priority, paper section, role/reason)
KEEP = {
 # --- Anchor / direct competitors ---
 '42749733': ('A', 'Intro; Methods; Results (all)', 'Study under test (Bahl et al. 2026)'),
 '37732899': ('A', 'Methods H1; Discussion', 'TNBC cell-line ATAC-seq (Harris 2023): candidate third lab for cross-lab H1'),
 '41539304': ('A', 'Methods H1; Discussion', 'Multi-omic breast cancer cell-line reference samples: reproducibility benchmark'),
 '42470318': ('A', 'Discussion (everolimus)', 'Competing mTOR-inhibitor sensitivity biomarker in TNBC (TOMM40/FADS2)'),
 # --- H3: TF dependencies (single-gene evidence for HC-TFs) ---
 '40640173': ('B', 'Results H3; Discussion', 'ID4 (HC-TF) as TNBC driver/target: prior dependency evidence'),
 '40586711': ('B', 'Results H3; Discussion', 'NFIB (HC-TF) regulation in TNBC'),
 '36171192': ('B', 'Results H3; Discussion', 'SNAI1 (HC-TF) loss phenotype in TNBC cells'),
 '29511313': ('B', 'Results H3; Discussion', 'FOXC1 (HC-TF) as therapeutic target in basal-like BC'),
 '36207293': ('B', 'Results H3; Discussion', 'MYC (HC-TF) regulatory axis in aggressive BC'),
 '39570746': ('B', 'Discussion', 'ZBTB7A/KDM5 TF networks in basal BC'),
 '29729076': ('C', 'Discussion', 'Snail1 genome-wide binding in TNBC cells'),
 '37217832': ('C', 'Discussion', 'Epigenetic reprogramming of TNBC cell states (dCas9)'),
 # --- TF activity -> drug response ---
 '42020604': ('A', 'Intro; Discussion', 'Master-regulator (VIPER-type) activity guiding combination therapy (Nat Genet 2026)'),
 '36374194': ('A', 'Intro', 'OncoLoop: network-based precision medicine framework (VIPER lineage)'),
 '39079653': ('B', 'Intro', 'TF activity for drug stratification (bladder)'),
 '37103222': ('B', 'Intro', 'Drug-target TF activity stratification (ovarian)'),
 '30704449': ('B', 'Intro', 'Pathway activity inference predicts drug sensitivity'),
 # --- Pharmacogenomic consistency (H4 rationale) ---
 '27905419': ('A', 'Intro; Discussion', 'Drug response consistency CCLE vs CGP (Nature 2016)'),
 '27193678': ('A', 'Intro; Discussion', 'Reproducible pharmacogenomic profiling (Nature 2016)'),
 '31846027': ('A', 'Intro; Methods H4', 'Evaluating consistency of large pharmacogenomic studies'),
 '30089904': ('A', 'Discussion', 'Cell-line genetic/transcriptional evolution alters drug response (Ben-David 2018)'),
 '29313488': ('A', 'Intro; Discussion', 'Replication study of GDSC drug-sensitivity markers (RP:CB)'),
 '35585485': ('A', 'Methods H4', 'Statistical approaches for association testing in noisy drug-screen data'),
 '29123200': ('B', 'Methods H4', 'Discordancy partitioning for inconsistent pharmacogenomic studies'),
 '27408686': ('B', 'Discussion', 'Assessment of pharmacogenomic agreement'),
 '28928933': ('B', 'Discussion', 'Revisiting inconsistency in pharmacogenomic studies'),
 '41726490': ('B', 'Discussion', 'Resolving inconsistency across pharmacogenomic datasets'),
 '41169498': ('B', 'Discussion', 'Spatial artefacts and drug-screen reproducibility'),
 '36259664': ('B', 'Methods H4', 'CREAMMIST integrative drug-response database'),
 '30553813': ('B', 'Methods H4', 'CellMinerCDB cross-database pharmacogenomics'),
 '37542042': ('B', 'Discussion', 'Proteomic profiling across breast cancer cell lines'),
 # --- Chromatin / ATAC methods (H1) ---
 '36452861': ('A', 'Methods H1', 'ATAC-seq normalisation benchmark (GC content)'),
 '32576878': ('B', 'Methods H1', 'Differential accessibility strategies compared'),
 '37990143': ('B', 'Methods H1', 'Reproducibility metrics for chromatin accessibility data'),
 '33878366': ('B', 'Methods H1', 'Reproducible peak sets from replicates (ChIP-R)'),
 '30791920': ('B', 'Methods H1', 'Protocol-specific bias in TF footprinting (ATAC vs DNase)'),
 '38817966': ('A', 'Methods H2', 'Benchmark of TF prioritisation tools'),
 '38429368': ('B', 'Discussion', 'FOXA1/GRHL2 from chromatin accessibility in breast cancer lines'),
 '29458327': ('C', 'Discussion', 'Histone modification profiling across breast cancer lines'),
 # --- Regulon / TF-activity methods (H2) ---
 '42741841': ('A', 'Methods H2; Discussion', 'TF regulatory atlas for activity inference (NAR 2026)'),
 '31340985': ('A', 'Methods H2', 'DoRothEA benchmark of TF-activity resources (Garcia-Alonso 2019)'),
 '42752548': ('B', 'Discussion', 'Benchmark of single-cell TF activity inference (2026)'),
 '32051003': ('B', 'Methods H2', 'Robustness of TF/pathway tools on scRNA-seq (Holland 2020)'),
 '39781510': ('B', 'Discussion', 'Context-specific TF regulons'),
 # --- CRISPR screen methods (H3) ---
 '39030569': ('A', 'Methods H3', 'Benchmark of bias correction in CRISPR-Cas9 screens'),
 # --- mTOR / everolimus context ---
 '27893038': ('A', 'Discussion (everolimus)', 'PI3K/AKT/mTOR targeting in mesenchymal TNBC (JAMA Oncol 2017)'),
 '40817986': ('A', 'Discussion (everolimus)', 'RCT carboplatin +/- everolimus, advanced TNBC'),
 '28270498': ('A', 'Discussion (everolimus)', 'RCT neoadjuvant cisplatin/paclitaxel +/- everolimus'),
 '24684785': ('B', 'Discussion (everolimus)', 'Phase 2 everolimus + carboplatin, metastatic TNBC'),
 '22703543': ('B', 'Discussion (everolimus)', 'Everolimus efficacy in basal-like TNBC cells'),
 '27374081': ('B', 'Discussion (everolimus)', 'mTOR targeting across TNBC molecular subtypes'),
 '28290700': ('B', 'Discussion (everolimus)', 'Metabolic response to everolimus in TNBC PDXs'),
 '24708766': ('B', 'Discussion (PDX)', 'TNBC PDXs reproduce tumour features and drug response'),
 '38211606': ('C', 'Discussion', 'FUTURE-SUPER subtype-based therapy trial (Lancet Oncol 2024)'),
 # --- Small-n PDX statistics (PDX re-test) ---
 '42970040': ('A', 'Methods (PDX re-test); Discussion', 'Multiple-comparison tests in small-sample PDX studies (2026)'),
 '38641411': ('A', 'Discussion (PDX)', 'NCI PDXNet consensus on PDX growth/response assessment'),
 '34107029': ('A', 'Discussion (PDX)', 'PDX design: number of lines vs subsampling'),
 '30944774': ('B', 'Discussion (PDX)', 'Factors influencing PDX response classification'),
 '33850213': ('B', 'Methods (PDX re-test)', 'Statistics for PDX tumour-growth experiments'),
 # --- Reproducibility framing ---
 '34874005': ('A', 'Intro', 'Investigating replicability of preclinical cancer biology (RP:CB)'),
 '34874008': ('A', 'Intro', 'Challenges for assessing replicability (RP:CB)'),
 '34874006': ('B', 'Intro', 'Is preclinical cancer research reproducible enough?'),
 '28100397': ('B', 'Discussion', 'Replication of drug-indication discovery from expression compendia'),
 '25490932': ('C', 'Intro', 'Launch of RP:CB'),
 '28662052': ('C', 'Discussion', 'Can researchers predict which reports will replicate?'),
 # --- TNBC cell-line models ---
 '21633166': ('A', 'Intro; Discussion', 'Lehmann 2011: TNBC subtypes and preclinical models'),
 '23151021': ('A', 'Methods; Discussion', 'Molecular characterisation of TNBC cell-line models'),
 '29348905': ('A', 'Discussion', 'Discrepancies between TNBC cell lines and tumours'),
 '28583138': ('B', 'Methods', 'Multi-omic characterisation of breast cancer cell lines'),
 '24162158': ('B', 'Methods', 'Cell lines and intrinsic subtypes'),
 '32353087': ('B', 'Discussion', 'Same TNBC lines classified differently by different methods'),
 '25344583': ('B', 'Discussion', 'Kinase dependency across TNBC subtypes'),
 '38391914': ('B', 'Discussion', 'Subtyping of breast cancer lines (Cells, MDPI 2024)'),
 '38892472': ('C', 'Intro', 'TNBC small-molecule therapeutics review (IJMS, MDPI 2024)'),
}

DROP_REASON = {
 'S3': 'Off-topic: "VIPER"/"protein activity"/"master regulator" used in unrelated sense or non-cancer',
 'S4': 'Prognostic-signature or cell-line-establishment paper; no reproducibility contribution',
 'S5': 'Single-cell integration or non-cancer ATAC method; not needed for bulk cross-lab comparison',
 'S6': 'Bacterial/plant regulatory networks or unrelated single-cell tool',
 'S7': 'Clinical/review paper without a predictive-biomarker angle',
 'S8a': 'Clinical prediction model using C-index; no small-sample methodology',
 'S8b': 'PDX application paper without design/statistics content',
 'S9': 'Individual replication study outside pharmacogenomics, or not about reproducibility',
 'S10': 'Single-compound cytotoxicity study in cell lines',
 'S1': 'Not breast cancer or not TF activity', 'S2': 'Single-gene mechanism unrelated to HC-TFs or not a TF',
}

rows = []
for _, r in df.iterrows():
    s = r.searches.split('|')
    k = KEEP.get(r.pmid)
    rows.append({
        'ID': f'PMID:{r.pmid}', 'Source': 'PubMed', 'Search': r.searches,
        'Citation': f"{r.first} et al. {r.yr}, {r.journal}" if int(r.nAu or 0) > 1 else f"{r.first} {r.yr}, {r.journal}",
        'Title': r.title, 'Year': r.yr, 'DOI': r.doi, 'PMID': r.pmid, 'PubType': r.pt,
        'Decision': 'keep' if k else 'drop', 'Priority': k[0] if k else '',
        'Reason / role': k[2] if k else DROP_REASON.get(s[0], 'Off-topic on title screen'),
        'Paper section': k[1] if k else '', 'Screen level': 'title (abstract if borderline)',
    })
log = pd.DataFrame(rows)
missing = set(KEEP) - set(log.PMID)
assert not missing, missing

# Citation chasing and methods references (verified in OpenAlex on 2026-11-05)
CHASE = [
 ('Bahl refs', 'Alvarez et al. 2016, Nat Genet', 'Functional characterization of somatic mutations in cancer using network-based inference of protein activity (VIPER)', 2016, '10.1038/ng.3593', 'A', 'Methods H2', 'VIPER method'),
 ('Bahl refs', 'Lachmann et al. 2016, Bioinformatics', 'ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information', 2016, '10.1093/bioinformatics/btw216', 'B', 'Methods H2', 'ARACNe network inference'),
 ('Bahl refs', 'Marcotte et al. 2016, Cell', 'Functional Genomic Landscape of Human Breast Cancer Drivers, Vulnerabilities, and Resistance', 2016, '10.1016/j.cell.2015.11.062', 'A', 'Methods H2/H3; Discussion', "Likely source of Bahl's cell-line RNA (GSE73526; to confirm) and earlier shRNA dependency screens in breast lines"),
 ('Bahl refs', 'Lehmann et al. 2021, Nat Commun', 'Multi-omics analysis identifies therapeutic vulnerabilities in triple-negative breast cancer subtypes', 2021, '10.1038/s41467-021-26502-6', 'A', 'Intro; Discussion', 'Closest prior work on TNBC subtype vulnerabilities'),
 ('Bahl refs', 'Corces et al. 2018, Science', 'The chromatin accessibility landscape of primary human cancers', 2018, '10.1126/science.aav1898', 'A', 'Intro; Methods H1', 'TCGA ATAC atlas (source of tumour ATAC)'),
 ('Bahl refs', 'Iorio et al. 2016, Cell', 'A Landscape of Pharmacogenomic Interactions in Cancer', 2016, '10.1016/j.cell.2016.06.017', 'A', 'Methods H4', 'GDSC landscape'),
 ('Bahl refs', 'Yang et al. 2013, Nucleic Acids Res', 'Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells', 2013, '10.1093/nar/gks1111', 'B', 'Methods H4', 'GDSC resource'),
 ('Bahl refs', 'Smirnov et al. 2016, Bioinformatics', 'PharmacoGx: an R package for analysis of large pharmacogenomic datasets', 2016, '10.1093/bioinformatics/btv723', 'A', 'Methods H4', "Bahl's pharmacogenomic toolkit (likely source of their C-index)"),
 ('Bahl refs', 'Gao et al. 2015, Nat Med', 'High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response', 2015, '10.1038/nm.3954', 'A', 'Discussion (PDX)', 'PDX screening design; single-animal response metrics'),
 ('Bahl refs', 'Lehmann et al. 2011, J Clin Invest', 'Identification of human triple-negative breast cancer subtypes and preclinical models for selection of targeted therapies', 2011, '10.1172/jci45014', 'A', 'Intro', 'TNBC subtypes (also PMID 21633166)'),
 ('Bahl refs', 'Seashore-Ludlow et al. 2015, Cancer Discov', 'Harnessing Connectivity in a Large-Scale Small-Molecule Sensitivity Dataset', 2015, '10.1158/2159-8290.cd-15-0235', 'B', 'Methods H4', 'CTRPv2'),
 ('Bahl refs', 'Mer et al. 2019, Cancer Res', 'Integrative Pharmacogenomics Analysis of Patient-Derived Xenografts', 2019, '10.1158/0008-5472.can-19-0349', 'A', 'Methods (PDX re-test)', 'Xeva/PDX pharmacogenomics: likely source of BestAvgResponse metric'),
 ('Bahl refs', 'Bahl et al. 2026, Code Ocean', 'Code capsule for Bahl et al. 2026', 2026, '10.24433/co.3288699.v1', 'A', 'Methods (PDX re-test)', 'Their code: inspect C-index implementation'),
 ('Bahl refs', 'Wang et al. 2022, iScience', 'NRF3 activates mTORC1 arginine-dependently for cancer cell viability', 2022, '10.1016/j.isci.2023.106045', 'B', 'Discussion (everolimus)', 'Biological rationale for NFE2L3-mTOR link'),
 ('Methods', 'Dempster et al. 2021, Genome Biol', 'Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects', 2021, '10.1186/s13059-021-02540-7', 'A', 'Methods H3', 'DepMap gene-effect scores'),
 ('Methods', 'Tsherniak et al. 2017, Cell', 'Defining a Cancer Dependency Map', 2017, '10.1016/j.cell.2017.06.010', 'A', 'Methods H3', 'DepMap'),
 ('Methods', 'Corsello et al. 2020, Nat Cancer', 'Discovering the anticancer potential of non-oncology drugs by systematic viability profiling', 2020, '10.1038/s43018-019-0018-6', 'A', 'Methods H4', 'PRISM Repurposing'),
 ('Methods', 'Haibe-Kains et al. 2013, Nature', 'Inconsistency in large pharmacogenomic studies', 2013, '10.1038/nature12831', 'A', 'Intro', 'Pharmacogenomic inconsistency'),
 ('Methods', 'Schep et al. 2017, Nat Methods', 'chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data', 2017, '10.1038/nmeth.4401', 'A', 'Methods H1', 'Motif deviations'),
 ('Methods', 'Li et al. 2011, Ann Appl Stat', 'Measuring reproducibility of high-throughput experiments', 2011, '10.1214/11-aoas466', 'A', 'Methods H1', 'IDR'),
 ('Methods', 'Harrell et al. 1996, Stat Med', 'Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors', 1996, '10.1002/(sici)1097-0258(19960229)15:4<361::aid-sim168>3.0.co;2-4', 'A', 'Methods (PDX re-test)', 'C-index definition'),
 ('Methods', 'Pencina & D\'Agostino 2004, Stat Med', 'Overall C as a measure of discrimination in survival analysis: model specific population value and confidence interval estimation', 2004, '10.1002/sim.1802', 'A', 'Methods (PDX re-test)', 'C-index CI estimation'),
 ('Methods', 'Schröder et al. 2011, Bioinformatics', 'survcomp: an R/Bioconductor package for performance assessment and comparison of survival models', 2011, '10.1093/bioinformatics/btr511', 'A', 'Methods (PDX re-test)', 'Likely C-index implementation behind PharmacoGx'),
 ('Methods', 'Phipson & Smyth 2010, Stat Appl Genet Mol Biol', 'Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn', 2010, '10.2202/1544-6115.1585', 'A', 'Methods; Discussion', 'Bahl report p = 0 for 58 drug-TF pairs'),
 ('Methods', 'Benjamini & Hochberg 1995, J R Stat Soc B', 'Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing', 1995, '10.1111/j.2517-6161.1995.tb02031.x', 'A', 'Methods', 'BH-FDR'),
 ('Methods', 'Badia-i-Mompel et al. 2022, Bioinform Adv', 'decoupleR: ensemble of computational methods to infer biological activities from omics data', 2022, '10.1093/bioadv/vbac016', 'B', 'Methods H2', 'Alternative activity-inference methods (sensitivity analysis)'),
]
chase = pd.DataFrame([{'ID': f'DOI:{c[4]}', 'Source': 'Citation chasing' if c[0] == 'Bahl refs' else 'Methods (hand-picked)',
    'Search': c[0], 'Citation': c[1], 'Title': c[2], 'Year': str(c[3]), 'DOI': c[4], 'PMID': '', 'PubType': '',
    'Decision': 'keep', 'Priority': c[5], 'Reason / role': c[7], 'Paper section': c[6],
    'Screen level': 'full record (OpenAlex)'} for c in CHASE])
log = pd.concat([log, chase], ignore_index=True)
log.to_csv(lit / 'literature_log_draft.csv', index=False)
k = log[log.Decision == 'keep']
print('records', len(log), '| keep', len(k), k.Priority.value_counts().to_dict(), '| drop', (log.Decision == 'drop').sum())
