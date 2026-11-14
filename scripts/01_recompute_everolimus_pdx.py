"""Re-test Bahl et al. 2026 (Nat Commun) everolimus PDX TF-activity associations
from their own Supplementary Data 8, using exact-permutation p-values (n = 11)."""
import itertools, numpy as np, pandas as pd, sys
from pathlib import Path
d = Path(sys.argv[1]); out = Path(sys.argv[2])
V = pd.read_csv(d/'SuppData8_VIPER_activity_PDXs.csv', index_col=0)
C = pd.read_csv(d/'SuppData8_everolimus_PDX_clinical.csv')
R = pd.read_csv(d/'SuppData8_everolimus_PDX_concordance.csv').set_index('TR')
ids = C['Sample.ID'].str.replace(r'T$', '', regex=True)
cols = [next(c for c in V.columns if c.startswith(i+'X')) for i in ids]
assert len(set(cols)) == 11
y = C['EVEROLIMUS_BestAvgResponse'].to_numpy(float)
I, J = np.triu_indices(len(y), 1)
def cidx(x, Y):  # Y: (B, n) permuted responses
    dx = np.sign(x[I] - x[J]); dy = np.sign(Y[:, I] - Y[:, J])
    s = dx * dy
    return ((s > 0) + 0.5 * (s == 0) * (dy != 0)).sum(1) / (dy != 0).sum(1)
# exact null: all 11! permutations is 39.9M; use 200k random permutations (seeded)
rng = np.random.default_rng(20261105)
P = np.array([rng.permutation(y) for _ in range(200_000)])
rows = []
for tf, r in V.iterrows():
    x = r[cols].to_numpy(float)
    c = cidx(x, y[None])[0]; null = np.abs(cidx(x, P) - 0.5)
    p = (1 + (null >= abs(c - 0.5) - 1e-12).sum()) / (1 + len(P))
    crep = R.loc[tf, 'c_index']
    p_at_rep = (1 + (null >= abs(crep - 0.5) - 1e-12).sum()) / (1 + len(P))
    rows.append(dict(TF=tf, c_recomputed=c, p_perm=p, c_reported=crep,
                     p_reported=R.loc[tf, 'p_value'], FDR_reported=R.loc[tf, 'FDR'],
                     p_perm_at_reported_c=p_at_rep))
df = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
df['FDR_perm'] = multipletests(df.p_perm, method='fdr_bh')[1]
df.sort_values('p_perm').to_csv(out, index=False, float_format='%.4g')
sig = df[df.FDR_reported <= 0.05].sort_values('p_reported')
print(sig[['TF','c_reported','p_reported','FDR_reported','c_recomputed','p_perm','p_perm_at_reported_c','FDR_perm']].to_string(index=False))
print('\nTFs reported FDR<=0.05:', len(sig), '| recomputed FDR_perm<=0.05:', (df.FDR_perm <= 0.05).sum(),
      '| min p_perm:', df.p_perm.min().round(4))
