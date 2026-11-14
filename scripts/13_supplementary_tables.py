"""Supplementary Tables S3-S7 as one .xlsx (values copied from results/; no recomputation)."""
import sys, pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font
root = sys.argv[1]; R = f'{root}/results'
sheets = {
 'S3_PDX_retest': [('Everolimus PDX re-test, 94 HC-TFs', pd.read_csv(f'{R}/01_everolimus_pdx_recomputed.csv'))],
 'S4_dependency': [('DepMap 23Q2 enrichment', pd.read_csv(f'{R}/05_h3_primary_h3a_enrichment.csv')),
                   ('Project Score enrichment', pd.read_csv(f'{R}/07ps_h3_primary_h3a_enrichment.csv')),
                   ('DepMap selectivity per TF (TNBC vs luminal)', pd.read_csv(f'{R}/05_h3a_selectivity_per_TF.csv')),
                   ('DepMap activity vs dependency per TF (H3b)', pd.read_csv(f'{R}/06_h3b_activity_dependency_per_TF.csv')),
                   ('Project Score activity vs dependency per TF', pd.read_csv(f'{R}/07ps_h3b_activity_dependency_per_TF.csv')),
                   ('Screen agreement (H3c, H3d)', pd.read_csv(f'{R}/08_h3c_h3d_agreement.csv')),
                   ('Subtypes (H3f)', pd.read_csv(f'{R}/09_h3f_subtypes.csv'))],
 'S5_PRISM': [('Summary (random-pair null from 200 permutations per pair; superseded by the row below)', pd.read_csv(f'{R}/10_h4_summary.csv')),
              ('Random-pair null, 10,000 permutations per pair (script 10b; reported in the main text)', pd.read_csv(f'{R}/10b_h4_null_10k_summary.csv')),
              ('All testable pairs', pd.read_csv(f'{R}/10_h4_prism_pairs.csv')),
              ('Pre-specified mTOR inhibitors x NFE2L3/CBFB', pd.read_csv(f'{R}/10_h4_mtor_prespecified.csv')),
              ('GDSC1 summary (script 16; second planned screen, same Sanger platform as GDSC1000)', pd.read_csv(f'{R}/16_h4_gdsc1_summary.csv')),
              ('GDSC1, all testable pairs', pd.read_csv(f'{R}/16_h4_gdsc1_pairs.csv')),
              ('GDSC1, mTOR-targeting compounds x NFE2L3/CBFB', pd.read_csv(f'{R}/16_h4_gdsc1_mtor_prespecified.csv'))],
 'S6_activity_reproducibility': [('Summary, DepMap 23Q2', pd.read_csv(f'{R}/11_h2_summary.csv')),
                                 ('Summary, Lund RNA-seq (GSE337499)', pd.read_csv(f'{R}/11b_h2_lund_summary.csv')),
                                 ('Per TF, DepMap', pd.read_csv(f'{R}/11_h2_depmap_activity_agreement.csv').rename(columns={'Unnamed: 0': 'TF'})),
                                 ('Per TF, Lund', pd.read_csv(f'{R}/11b_h2_lund_activity_agreement.csv').rename(columns={'Unnamed: 0': 'TF'})),
                                 ('Line identity, Lund (exploratory; rows released, columns Lund)', pd.read_csv(f'{R}/11b_h2_lund_line_identity.csv'))],
 'S7_crosslab_chromatin': [('Criterion and identity tests per shared line', pd.read_csv(f'{R}/14_h1_crosslab.csv')),
                           ('Sensitivity analysis: called/not-called regions, Jaccard index (script 14b)', pd.read_csv(f'{R}/14b_h1_binarised.csv')),
                           ('Replicate agreement', pd.read_csv(f'{R}/14_h1_replicates.csv')),
                           ('Shared lines, raw (main-text Figure 2A)', pd.read_csv(f'{R}/14_h1_shared_lines_corr_raw.csv')),
                           ('Shared lines, centred (main-text Figure 2B)', pd.read_csv(f'{R}/14_h1_shared_lines_corr_centred.csv'))],
}
out = f'{root}/supplementary/Supplementary_Tables_S3-S7.xlsx'; heads = []
with pd.ExcelWriter(out, engine='openpyxl') as w:
    for sh, blocks in sheets.items():
        r = 0
        for title, df in blocks:
            pd.DataFrame([[title]]).to_excel(w, sheet_name=sh, startrow=r, index=False, header=False); heads.append((sh, r + 1))
            df.to_excel(w, sheet_name=sh, startrow=r + 1, index=False); r += len(df) + 4
wb = load_workbook(out)
for ws in wb:
    for row in ws.iter_rows():
        for c in row: c.font = Font(name='Arial', size=10)
for sh, r in heads: wb[sh].cell(r, 1).font = Font(name='Arial', size=11, bold=True)
wb.save(out); print('written', out)
