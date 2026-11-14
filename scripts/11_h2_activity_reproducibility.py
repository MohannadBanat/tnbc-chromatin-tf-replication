"""H2 (plan 2026-11-05): reproducibility of HC-TF activity in independent expression data.

Reference: Bahl et al.'s released cell-line activities (Supp. Data 5; computed on GSE73526 RNA-seq).
Test: VIPER (decoupler) with their full TCGA regulons on an independent expression source for the same
lines (here DepMap 23Q2 RNA-seq; Lund GSE337499 when available), genes z-scored across the shared lines.
Per TF: Spearman across shared lines between released and recomputed activity; median over HC-TFs with
bootstrap 95% CI (resampling lines). Null: line-label permutation (1,000), because released activities
exist only for 150 TNBC-specific TFs, so the plan's size-matched-regulon null cannot be formed
(logged as a deviation).
"""
import sys, re, numpy as np, pandas as pd, decoupler as dc
from scipy.stats import spearmanr
root = sys.argv[1]; R = f'{root}/data/depmap_raw'
rng = np.random.default_rng(20261106)
norm = lambda s: re.sub(r'[^A-Z0-9]', '', str(s).upper())
rel = pd.read_csv(f'{root}/data/bahl2026/SuppData5_TNBC_TF_activity_cell_lines.csv')
rel = rel.drop(columns=[rel.columns[0]]).set_index(rel.columns[1])   # rows TFs, cols lines
hc = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())
net = pd.read_csv(f'{root}/data/bahl2026/SuppData2_TCGA_VIPER_network_full.csv'); net.columns = ['source', 'target', 'weight']
net = net.groupby(['source', 'target'], as_index=False).weight.mean()

m = pd.read_csv(f'{R}/depmap23q2_Model.csv')
name2id = {norm(n): i for n, i in zip(m.StrippedCellLineName, m.ModelID)}
ex = pd.read_csv(f'{R}/depmap23q2_OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
ex.columns = [re.sub(r' \(\d+\)$', '', c) for c in ex.columns]; ex = ex.T.groupby(level=0).mean().T
shared = {c: name2id[norm(c)] for c in rel.columns if norm(c) in name2id and name2id[norm(c)] in ex.index}
print(f'Bahl cell lines: {rel.shape[1]}; matched to DepMap expression: {len(shared)}')
cols = list(shared); ids = [shared[c] for c in cols]
X = ex.loc[ids]; X = X.loc[:, X.std() > 0]; Z = (X - X.mean()) / X.std()
res = dc.mt.viper(Z, net[net.target.isin(Z.columns)], tmin=5); nes = res[0] if isinstance(res, tuple) else res
nes.index = cols
tfs = [t for t in rel.index if t in nes.columns]
A = rel.loc[tfs, cols].T.astype(float); Bm = nes[tfs]
rho = pd.Series({t: spearmanr(A[t], Bm[t]).statistic for t in tfs})
is_hc = rho.index.isin(hc)
def med_hc(Aa, Bb):
    return np.nanmedian([spearmanr(Aa[t], Bb[t]).statistic for t in tfs if t in hc])
obs = rho[is_hc].median()
boot = []
for _ in range(1000):
    idx = rng.choice(len(cols), len(cols)); boot.append(med_hc(A.iloc[idx], Bm.iloc[idx]))
perm = []
for _ in range(1000):
    pi = rng.permutation(len(cols)); perm.append(med_hc(A.iloc[pi].set_axis(A.index), Bm))
p = (1 + np.sum(np.array(perm) >= obs)) / 1001
print(f'HC-TFs compared: {is_hc.sum()}; median per-TF rho = {obs:.3f} (bootstrap 95% CI {np.nanpercentile(boot,2.5):.3f}-{np.nanpercentile(boot,97.5):.3f}); '
      f'permutation null median {np.median(perm):.3f} (95% {np.percentile(perm,2.5):.3f}-{np.percentile(perm,97.5):.3f}); one-sided p = {p:.4f}')
print(f'HC-TFs with rho >= 0.5: {(rho[is_hc] >= 0.5).sum()}; rho < 0.2: {(rho[is_hc] < 0.2).sum()}; non-HC TNBC-specific TFs median rho = {rho[~is_hc].median():.3f}')
out = pd.DataFrame({'rho_released_vs_DepMap': rho, 'HC': is_hc}); out.to_csv(f'{root}/results/11_h2_depmap_activity_agreement.csv', float_format='%.4g')
pd.DataFrame([dict(source='DepMap 23Q2', n_lines=len(cols), n_HC=int(is_hc.sum()), median_rho=obs,
                   ci_lo=np.nanpercentile(boot, 2.5), ci_hi=np.nanpercentile(boot, 97.5), perm_median=np.median(perm), p=p)]
             ).to_csv(f'{root}/results/11_h2_summary.csv', index=False)
for t in ['NFE2L3', 'CBFB', 'MYC', 'SOX9', 'FOXC1', 'E2F3']:
    if t in rho: print(f'  {t}: rho = {rho[t]:.3f}')
