"""Freeze cell-line groups from DepMap 23Q2 metadata ONLY (Model.csv + CRISPRScreenMap.csv),
before any gene-effect, expression or copy-number value is read (plan amendment 01, H3/H3a/H3f).

Rules (fixed here, applied mechanically):
- Breast line: CCLEName ends with '_BREAST'.
- TNBC: breast line with LegacySubSubtype == 'ERneg_HER2neg' (receptor-based, clinical definition).
- Luminal comparator: breast line with LegacySubSubtype == 'ERpos_HER2neg'.
- Other lineages: all non-breast lines (for the pan-cancer selectivity contrast).
- TNBC subtypes (H3f, descriptive), from LegacyMolecularSubtype:
    basal_A  -> 'basal_A';  basal_B -> 'basal_B';  'basal' (unspecified) -> 'basal_unspecified';
    luminal / HER2_amp / luminal_HER2_amp within TNBC -> 'luminal_like_TNBC' (proxy for LAR).
- Analysed lines = lines with a CRISPR (Avana/KY) screen in CRISPRScreenMap.
"""
import sys, hashlib, pandas as pd
d = sys.argv[1]; out = sys.argv[2]
m = pd.read_csv(f'{d}/depmap23q2_Model.csv'); m.columns = [c.strip() for c in m.columns]
screened = set(pd.read_csv(f'{d}/depmap23q2_CRISPRScreenMap.csv').ModelID)
m['is_breast'] = m.CCLEName.fillna('').str.endswith('_BREAST')
def group(r):
    if not r.is_breast: return 'other_lineage'
    if r.LegacySubSubtype == 'ERneg_HER2neg': return 'TNBC'
    if r.LegacySubSubtype == 'ERpos_HER2neg': return 'luminal'
    return 'breast_other'
m['group'] = m.apply(group, axis=1)
sub = {'basal_A': 'basal_A', 'basal_B': 'basal_B', 'basal': 'basal_unspecified',
       'luminal': 'luminal_like_TNBC', 'HER2_amp': 'luminal_like_TNBC', 'luminal_HER2_amp': 'luminal_like_TNBC'}
m['tnbc_subtype'] = [sub.get(s, 'unannotated') if g == 'TNBC' else '' for s, g in zip(m.LegacyMolecularSubtype, m.group)]
m['crispr_screened'] = m.ModelID.isin(screened)
cols = ['ModelID', 'StrippedCellLineName', 'group', 'tnbc_subtype', 'LegacyMolecularSubtype', 'LegacySubSubtype', 'crispr_screened']
m[cols].sort_values(['group', 'tnbc_subtype', 'StrippedCellLineName']).to_csv(out, index=False)
s = m[m.crispr_screened]
print('screened lines by group:', s.group.value_counts().to_dict())
print('screened TNBC by subtype:', s[s.group == 'TNBC'].tnbc_subtype.value_counts().to_dict())
print('TNBC screened:', ', '.join(sorted(s[s.group == 'TNBC'].StrippedCellLineName)))
print('luminal screened:', ', '.join(sorted(s[s.group == 'luminal'].StrippedCellLineName)))
print('sha256', hashlib.sha256(open(out, 'rb').read()).hexdigest())
