"""H3e confounders and H3f subtype structure (amendment 01), DepMap 23Q2.
H3e: no TNBC-selective dependency was found (script 05), so the per-dependency copy-number test has no
targets. Applied instead to the only FDR-significant activity-dependency association from H3b (ZIC4):
partial Spearman adjusting for growth rate, screen quality and the TF's copy number.
H3f: per TNBC subtype (frozen), HC-TFs with median gene effect <= -0.5, split by common-essential status.
"""
import sys, re, numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
root = sys.argv[1]; R = f'{root}/data/depmap_raw'
g = pd.read_csv(f'{root}/groups/cell_line_groups_FROZEN_2026-11-06.csv')
hc = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())
hdr = pd.read_csv(f'{R}/depmap23q2_CRISPRGeneEffect.csv', nrows=0).columns
sym = {c: re.sub(r' \(\d+\)$', '', c) for c in hdr[1:]}
ge = pd.read_csv(f'{R}/depmap23q2_CRISPRGeneEffect.csv', index_col=0, usecols=[hdr[0]] + [c for c in hdr[1:] if sym[c] in set(hc)])
ge.columns = [sym[c] for c in ge.columns]
ce = set(pd.read_csv(f'{R}/depmap23q2_CRISPRInferredCommonEssentials.csv').iloc[:, 0].str.replace(r' \(\d+\)$', '', regex=True))

# H3f
tn = g[(g.group == 'TNBC') & g.ModelID.isin(ge.index)]
rows = []
for st, sub in tn.groupby('tnbc_subtype'):
    med = ge.loc[sub.ModelID].median()
    dep = med[med <= -0.5]
    rows.append(dict(subtype=st, n_lines=len(sub), n_dep=len(dep),
                     non_common_essential=', '.join(sorted(t for t in dep.index if t not in ce)) or '-',
                     common_essential=', '.join(sorted(t for t in dep.index if t in ce))))
h3f = pd.DataFrame(rows); h3f.to_csv(f'{root}/results/09_h3f_subtypes.csv', index=False)
print('H3f subtypes:'); print(h3f.to_string(index=False))

# H3e on ZIC4
b = pd.read_csv(f'{root}/results/06_h3b_activity_dependency_per_TF.csv')
hits = b[(b.HC) & (b.FDR_activity_HC <= 0.05)].TF.tolist()
print('\nH3b FDR-significant HC-TFs:', hits)
gr = pd.read_csv(f'{R}/depmap23q2_CRISPRInferredModelGrowthRate.csv', index_col=0).iloc[:, 0]
qc = pd.read_csv(f'{R}/depmap23q2_AchillesScreenQCReport.csv')
smap = pd.read_csv(f'{R}/depmap23q2_CRISPRScreenMap.csv')
qcol = [c for c in qc.columns if 'NNMD' in c.upper() or 'nnmd' in c] or [c for c in qc.columns if 'ROC' in c.upper()]
q = qc.groupby('ModelID')['ScreenNNMD'].mean(); qcol = ['ScreenNNMD']
cnh = pd.read_csv(f'{R}/depmap23q2_OmicsCNGene.csv', nrows=0).columns
cn = pd.read_csv(f'{R}/depmap23q2_OmicsCNGene.csv', index_col=0, usecols=[cnh[0]] + [c for c in cnh[1:] if re.sub(r' \(\d+\)$', '', c) in set(hits)])
cn.columns = [re.sub(r' \(\d+\)$', '', c) for c in cn.columns]
import decoupler as dc  # noqa: F401  (activity re-read from results not stored per line; recompute minimal)
for t in hits:
    print(f'{t}: H3b rho = {b.loc[b.TF == t, "rho_activity"].iloc[0]:.3f}; median gene effect in TNBC = {b.loc[b.TF == t, "median_effect"].iloc[0]:.3f}')
    lines = [m for m in tn.ModelID if m in gr.index and m in q.index and m in cn.index]
    cov = pd.DataFrame({'growth': gr[lines], 'qc': q[lines], 'cn': cn.loc[lines, t]})
    cov = cov.dropna(); lines = list(cov.index); eff = ge.loc[lines, t]
    # technical-only model: rank-regress gene effect on covariates
    Xr = np.column_stack([np.ones(len(lines))] + [rankdata(cov[c]) for c in cov])
    beta, *_ = np.linalg.lstsq(Xr, rankdata(eff), rcond=None)
    r2 = 1 - np.sum((rankdata(eff) - Xr @ beta) ** 2) / np.sum((rankdata(eff) - rankdata(eff).mean()) ** 2)
    print(f'  technical-only model (growth, screen QC [{qcol[0] if qcol else "n/a"}], copy number) explains R2 = {r2:.2f} of the ranked gene effect (n = {len(lines)})')
    print(f'  copy number vs gene effect: rho = {spearmanr(cov.cn, eff).statistic:.3f}')
