"""H3b (amendment 01): does HC-TF activity predict dependency within TNBC lines?

- Activity: VIPER (decoupler 2.x, dc.mt.viper) with Bahl et al.'s released TCGA network (all 1,996
  regulons; signed VIPER weights as mode of regulation), on DepMap 23Q2 expression of the frozen TNBC
  lines, each gene z-scored across those lines.
- Per TF: Spearman rho between activity (NES) and Chronos gene effect across TNBC lines
  (negative rho = more active lines depend more on the TF).
- HC-TF summary: median rho; one-sample sign test of rho < 0; comparison with 1,000 regulon-size-matched
  draws of non-HC TFs (two-sided empirical p).
- Secondary: same with the TF's own mRNA; paired Wilcoxon of (rho_activity - rho_mRNA).
"""
import sys, re, numpy as np, pandas as pd, decoupler as dc
from scipy.stats import spearmanr, binomtest, wilcoxon
from statsmodels.stats.multitest import multipletests
root = sys.argv[1]; R = f'{root}/data/depmap_raw'
GE = sys.argv[2] if len(sys.argv) > 2 else 'depmap23q2_CRISPRGeneEffect.csv'
TAG = sys.argv[3] if len(sys.argv) > 3 else '06'
rng = np.random.default_rng(20261106)

g = pd.read_csv(f'{root}/groups/cell_line_groups_FROZEN_2026-11-06.csv')
net = pd.read_csv(f'{root}/data/bahl2026/SuppData2_TCGA_VIPER_network_full.csv')
net.columns = ['source', 'target', 'weight']
net = net.groupby(['source', 'target'], as_index=False).weight.mean()
hc = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())

ex = pd.read_csv(f'{R}/depmap23q2_OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
ex.columns = [re.sub(r' \(\d+\)$', '', c) for c in ex.columns]
ex = ex.T.groupby(level=0).mean().T
gehdr = pd.read_csv(f'{R}/{GE}', nrows=0).columns
sym = {c: re.sub(r' \(\d+\)$', '', c) for c in gehdr[1:]}
tfs = set(net.source)
ge = pd.read_csv(f'{R}/{GE}', index_col=0,
                 usecols=[gehdr[0]] + [c for c in gehdr[1:] if sym[c] in tfs])
ge.columns = [sym[c] for c in ge.columns]; ge = ge.T.groupby(level=0).mean().T

tn = [m for m in g.loc[g.group == 'TNBC', 'ModelID'] if m in ex.index and m in ge.index]
X = ex.loc[tn]
X = X.loc[:, X.std() > 0]
Z = (X - X.mean()) / X.std()
res = dc.mt.viper(Z, net[net.target.isin(Z.columns)], tmin=5)
nes = res[0] if isinstance(res, tuple) else res
print(f'TNBC lines: {len(tn)}; genes in signature: {Z.shape[1]}; TFs with activity: {nes.shape[1]}')

rows = []
for t in nes.columns:
    if t not in ge.columns or ge.loc[tn, t].isna().any():
        continue
    dep = ge.loc[tn, t]
    ra, pa = spearmanr(nes[t], dep)
    rm, pm = spearmanr(X[t], dep) if t in X.columns else (np.nan, np.nan)
    rows.append(dict(TF=t, HC=t in hc, n_targets=int((net.source == t).sum()),
                     rho_activity=ra, p_activity=pa, rho_mRNA=rm, p_mRNA=pm,
                     median_effect=dep.median(), sd_effect=dep.std()))
d = pd.DataFrame(rows)
d['FDR_activity_HC'] = np.nan
m = d.HC
d.loc[m, 'FDR_activity_HC'] = multipletests(d.loc[m, 'p_activity'], method='fdr_bh')[1]
d.to_csv(f'{root}/results/{TAG}_h3b_activity_dependency_per_TF.csv', index=False, float_format='%.4g')

h = d[d.HC].dropna(subset=['rho_activity'])
obs = h.rho_activity.median()
sign = binomtest(int((h.rho_activity < 0).sum()), len(h), 0.5, alternative='greater')
# size-matched null from non-HC TFs
pool = d[~d.HC].dropna(subset=['rho_activity']).copy()
bins = np.unique(np.quantile(d.n_targets, np.linspace(0, 1, 11)))
pool['bin'] = np.digitize(pool.n_targets, bins[1:-1]); hb = np.digitize(h.n_targets, bins[1:-1])
need = pd.Series(hb).value_counts()
byb = {b: pool.loc[pool.bin == b, 'rho_activity'].to_numpy() for b in need.index}
null = np.array([np.median(np.concatenate([rng.choice(byb[b], n, replace=len(byb[b]) < n) for b, n in need.items()]))
                 for _ in range(1000)])
e = null.mean(); p_null = (1 + np.sum(np.abs(null - e) >= abs(obs - e) - 1e-12)) / 1001
diff = (h.rho_activity - h.rho_mRNA).dropna()
w = wilcoxon(diff)
summary = dict(n_HC=len(h), median_rho_activity=obs, n_negative=int((h.rho_activity < 0).sum()),
               sign_test_p_one_sided=sign.pvalue, null_median=e, null_lo=np.percentile(null, 2.5),
               null_hi=np.percentile(null, 97.5), p_vs_null=p_null,
               median_rho_mRNA=h.rho_mRNA.median(), median_paired_diff=diff.median(), wilcoxon_p=w.pvalue,
               n_HC_FDR05=int((h.FDR_activity_HC <= 0.05).sum()))
pd.DataFrame([summary]).to_csv(f'{root}/results/{TAG}_h3b_summary.csv', index=False)
for k, v in summary.items(): print(f'{k}: {v:.4g}' if isinstance(v, float) else f'{k}: {v}')
pd.set_option('display.width', 200)
print('\nHC-TFs, most negative rho (activity):')
print(h.sort_values('rho_activity').head(10)[['TF', 'rho_activity', 'p_activity', 'FDR_activity_HC', 'rho_mRNA', 'median_effect']].to_string(index=False, float_format='%.3g'))
print('\nKey TFs:'); print(h[h.TF.isin(['NFE2L3', 'CBFB', 'MYC', 'SOX9', 'E2F3', 'FOXC1', 'ID4'])][['TF', 'rho_activity', 'p_activity', 'rho_mRNA', 'median_effect']].to_string(index=False, float_format='%.3g'))
