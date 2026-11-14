"""H3 primary + H3a (plan 2026-11-05 + amendment 01). DepMap 23Q2 Chronos gene effect.

Primary: proportion of HC-TFs with median gene effect <= -0.5 across frozen TNBC lines, vs
10,000 expression-decile-matched draws of non-HC TFs from the 1,996-TF TCGA network.
H3a: (i) same test excluding DepMap common essentials (CRISPRInferredCommonEssentials);
(ii) per-HC-TF selectivity: TNBC vs luminal, TNBC vs other lineages (Mann-Whitney, BH-FDR),
difference in medians with 2,000-sample bootstrap 95% CI; (iii) TNBC-selective dependency =
median TNBC <= -0.5 AND not common-essential AND TNBC<luminal at FDR <= 0.05.
"""
import sys, re, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
root = sys.argv[1]
GE = sys.argv[2] if len(sys.argv) > 2 else 'depmap23q2_CRISPRGeneEffect.csv'
CE = sys.argv[3] if len(sys.argv) > 3 else 'depmap23q2_CRISPRInferredCommonEssentials.csv'
TAG = sys.argv[4] if len(sys.argv) > 4 else '05'
R = f'{root}/data/depmap_raw'
rng = np.random.default_rng(20261106)
B = 10_000

g = pd.read_csv(f'{root}/groups/cell_line_groups_FROZEN_2026-11-06.csv')
hc = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())
universe = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData2_TCGA_VIPER_network_full.csv').iloc[:, 0].unique())

def load(path, genes):
    hdr = pd.read_csv(path, nrows=0).columns
    sym = {c: re.sub(r' \(\d+\)$', '', c) for c in hdr[1:]}
    keep = [hdr[0]] + [c for c in hdr[1:] if sym[c] in genes]
    df = pd.read_csv(path, usecols=keep, index_col=0)
    df.columns = [sym[c] for c in df.columns]
    return df.T.groupby(level=0).mean().T   # collapse any duplicated symbols

ge = load(f'{R}/{GE}', set(universe))
ex = load(f'{R}/depmap23q2_OmicsExpressionProteinCodingGenesTPMLogp1.csv', set(universe))
ce_raw = pd.read_csv(f'{R}/{CE}').iloc[:, 0]
common_ess = set(ce_raw.str.replace(r' \(\d+\)$', '', regex=True))

ids = lambda grp: [m for m in g.loc[g.group == grp, 'ModelID'] if m in ge.index]
tn, lu, ot = ids('TNBC'), ids('luminal'), ids('other_lineage')
tn_ex = [m for m in tn if m in ex.index]
genes = [t for t in universe if t in ge.columns and t in ex.columns]
hc_in = [t for t in hc if t in genes]
print(f'lines with gene effect: TNBC {len(tn)}, luminal {len(lu)}, other {len(ot)}; TNBC with expression {len(tn_ex)}')
print(f'TFs testable: {len(genes)} of {len(universe)}; HC-TFs testable: {len(hc_in)} of {len(hc)}; missing: {sorted(set(hc)-set(hc_in))}')

med = ge.loc[tn, genes].median()
expr = ex.loc[tn_ex, genes].median()
dec = pd.qcut(expr.rank(method='first'), 10, labels=False)

def matched_test(hc_set, pool, label):
    obs = (med[hc_set] <= -0.5).mean()
    pool = [t for t in pool if t not in hc_set]
    by_dec = {d: np.array([t for t in pool if dec[t] == d]) for d in range(10)}
    need = pd.Series(dec[hc_set]).value_counts()
    medv = med.to_dict(); null = np.empty(B)
    for b in range(B):
        draw = np.concatenate([rng.choice(by_dec[d], size=n, replace=len(by_dec[d]) < n) for d, n in need.items()])
        null[b] = np.mean([medv[t] <= -0.5 for t in draw])
    e = null.mean()
    p = (1 + np.sum(np.abs(null - e) >= abs(obs - e) - 1e-12)) / (1 + B)
    lo, hi = np.percentile(null, [2.5, 97.5])
    print(f'{label}: HC-TF dependency rate {obs:.3f} ({int(round(obs*len(hc_set)))}/{len(hc_set)}) | '
          f'matched null mean {e:.3f} [95% {lo:.3f}-{hi:.3f}] | fold {obs/e if e else float("nan"):.2f} | two-sided p {p:.4f}')
    return dict(analysis=label, n_hc=len(hc_set), n_dep=int(round(obs*len(hc_set))), rate=obs, null_mean=e,
                null_lo=lo, null_hi=hi, fold=obs/e if e else np.nan, p=p)

res = [matched_test(hc_in, genes, 'H3 primary (all TFs)')]
hc_nce = [t for t in hc_in if t not in common_ess]
res.append(matched_test(hc_nce, [t for t in genes if t not in common_ess], 'H3a excluding common essentials'))
pd.DataFrame(res).to_csv(f'{root}/results/{TAG}_h3_primary_h3a_enrichment.csv', index=False)

# H3a selectivity per HC-TF
def boot_diff(a, b, n=2000):
    a, b = np.asarray(a), np.asarray(b)
    d = [np.median(rng.choice(a, a.size)) - np.median(rng.choice(b, b.size)) for _ in range(n)]
    return np.percentile(d, [2.5, 97.5])
rows = []
for t in hc_in:
    a, l, o = ge.loc[tn, t].dropna(), ge.loc[lu, t].dropna(), ge.loc[ot, t].dropna()
    p_l = mannwhitneyu(a, l, alternative='two-sided').pvalue
    p_o = mannwhitneyu(a, o, alternative='two-sided').pvalue
    cl, cu = boot_diff(a, l)
    rows.append(dict(TF=t, median_TNBC=a.median(), median_luminal=l.median(), median_other=o.median(),
                     diff_TNBC_minus_luminal=a.median() - l.median(), diff_CI_lo=cl, diff_CI_hi=cu,
                     p_vs_luminal=p_l, p_vs_other=p_o, common_essential=t in common_ess,
                     median_expr_TNBC=expr[t]))
s = pd.DataFrame(rows)
s['FDR_vs_luminal'] = multipletests(s.p_vs_luminal, method='fdr_bh')[1]
s['FDR_vs_other'] = multipletests(s.p_vs_other, method='fdr_bh')[1]
s['dependency_TNBC'] = s.median_TNBC <= -0.5
s['TNBC_selective'] = s.dependency_TNBC & ~s.common_essential & (s.diff_TNBC_minus_luminal < 0) & (s.FDR_vs_luminal <= 0.05)
s = s.sort_values('median_TNBC')
s.to_csv(f'{root}/results/{TAG}_h3a_selectivity_per_TF.csv', index=False, float_format='%.4g')
pd.set_option('display.width', 220)
print('\nHC-TFs that are TNBC dependencies (median <= -0.5):')
print(s[s.dependency_TNBC][['TF', 'median_TNBC', 'median_luminal', 'median_other', 'common_essential', 'FDR_vs_luminal', 'FDR_vs_other', 'TNBC_selective']].to_string(index=False, float_format='%.3g'))
print('\nTNBC-selective dependencies (all three criteria):', s[s.TNBC_selective].TF.tolist())
print('Nominally TNBC < luminal (p<0.05, any dependency level):', s[(s.p_vs_luminal < 0.05) & (s.diff_TNBC_minus_luminal < 0)].TF.tolist())
