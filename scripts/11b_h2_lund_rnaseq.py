"""H2 part (i): Lund RNA-seq (GSE337499, FPKM, 8 TNBC lines x 2 replicates) vs Bahl et al.'s released activities.
Same method as script 11: log2(FPKM+1), replicate mean, genes collapsed by symbol, z-scored across the shared
lines, VIPER (decoupler, tmin=5) with Bahl's full TCGA regulons. Per-TF Spearman across shared lines; median
over HC-TFs; null = all distinct line permutations (exact when n! is small). Secondary (exploratory):
line identity - per-line Spearman across all released TFs, same line vs different lines.
"""
import sys, re, itertools, math, numpy as np, pandas as pd, decoupler as dc
from scipy.stats import spearmanr
root = sys.argv[1]
norm = lambda s: re.sub(r'[^A-Z0-9]', '', str(s).upper())
rel = pd.read_csv(f'{root}/data/bahl2026/SuppData5_TNBC_TF_activity_cell_lines.csv')
rel = rel.drop(columns=[rel.columns[0]]).set_index(rel.columns[1])
hc = set(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())
net = pd.read_csv(f'{root}/data/bahl2026/SuppData2_TCGA_VIPER_network_full.csv'); net.columns = ['source', 'target', 'weight']
net = net.groupby(['source', 'target'], as_index=False).weight.mean()
f = pd.read_csv(f'{root}/data/geo/GSE337499/GSE337499_fpkm_matrix.tsv.gz', sep='\t')
f = f[f.BT549_rep1 != 'FPKM'].copy()   # one repeated header row inside the file
num = [c for c in f.columns if c not in ('gene_id', 'gene_name', 'transcript_ids')]; f[num] = f[num].astype(float)
v = np.log2(f.drop(columns=['gene_id', 'gene_name', 'transcript_ids']) + 1)
v.columns = [re.sub(r'_rep\d+$', '', c) for c in v.columns]
v = v.T.groupby(level=0).mean().T; v.index = f.gene_name; v = v.groupby(level=0).mean()
relmap = {norm(c): c for c in rel.columns}
shared = [c for c in v.columns if norm(c) in relmap]
print(f'Lund lines: {v.shape[1]}; shared with released activities: {len(shared)} ({", ".join(shared)})')
X = v[shared].T; X = X.loc[:, X.std() > 0]; Z = (X - X.mean()) / X.std()
res = dc.mt.viper(Z, net[net.target.isin(Z.columns)], tmin=5); nes = res[0] if isinstance(res, tuple) else res
tfs = [t for t in rel.index if t in nes.columns]
A = rel.loc[tfs, [relmap[norm(c)] for c in shared]].T.astype(float); A.index = shared; B = nes[tfs]
rho = pd.Series({t: spearmanr(A[t], B[t]).statistic for t in tfs}); is_hc = rho.index.isin(hc)
obs = np.nanmedian(rho[is_hc]); n = len(shared)
perms = [p for p in itertools.permutations(range(n)) if p != tuple(range(n))]
null = np.array([np.nanmedian([spearmanr(A[t].to_numpy()[list(p)], B[t]).statistic for t in tfs if t in hc]) for p in perms])
pval = (1 + np.sum(null >= obs)) / (len(null) + 1)
rng = np.random.default_rng(20261106); boot = []
for _ in range(1000):
    idx = rng.choice(n, n)
    if len(set(idx)) < 3: continue
    boot.append(np.nanmedian([spearmanr(A[t].iloc[idx], B[t].iloc[idx]).statistic for t in tfs if t in hc]))
print(f'HC-TFs: {is_hc.sum()}; median per-TF rho = {obs:.3f} (bootstrap 95% CI {np.nanpercentile(boot,2.5):.3f}-{np.nanpercentile(boot,97.5):.3f}); '
      f'exact permutation null ({len(null)} permutations) median {np.median(null):.3f}; one-sided p = {pval:.4f}')
print(f'HC-TFs rho >= 0.5: {(rho[is_hc] >= 0.5).sum()}; non-HC median {np.nanmedian(rho[~is_hc]):.3f}')
for t in ['NFE2L3', 'CBFB']: print(f'  {t}: rho = {rho.get(t, np.nan):.3f}')
# exploratory line identity
M = pd.DataFrame({(a, b): [spearmanr(A.loc[a], B.loc[b]).statistic] for a in shared for b in shared}).T[0].unstack()
same = np.diag(M.to_numpy()); top = (M.to_numpy().argmax(1) == np.arange(n)).sum()
print(f'Line identity (exploratory): same-line rho median {np.median(same):.3f}, different-line median {np.median(M.to_numpy()[~np.eye(n, dtype=bool)]):.3f}; '
      f'released profile best matched by its own Lund line for {top}/{n} lines')
pd.DataFrame({'rho_released_vs_Lund': rho, 'HC': is_hc}).to_csv(f'{root}/results/11b_h2_lund_activity_agreement.csv', float_format='%.4g')
M.to_csv(f'{root}/results/11b_h2_lund_line_identity.csv', float_format='%.4g')
pd.DataFrame([dict(source='Lund GSE337499', n_lines=n, n_HC=int(is_hc.sum()), median_rho=obs, ci_lo=np.nanpercentile(boot, 2.5),
                   ci_hi=np.nanpercentile(boot, 97.5), perm_median=np.median(null), p=pval, n_perm=len(null), line_identity_top1=top)]
             ).to_csv(f'{root}/results/11b_h2_lund_summary.csv', index=False)
