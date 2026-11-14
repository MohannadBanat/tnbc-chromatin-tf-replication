"""H4 (plan 2026-11-05): replication of Bahl et al.'s cross-dataset drug-TF pairs in PRISM Repurposing 20Q2
(secondary screen, dose-response AUC).

- Pairs: drug-TF pairs significant (their FDR <= 0.05) in >= 2 of GRAY/GDSC2/CTRPv2 (Supp. Data 7).
  Reported direction = sign of (mean reported C - 0.5).
- Lines: frozen TNBC group with PRISM AUC and DepMap 23Q2 expression. Sensitivity = 1 - AUC
  (Bahl et al. use AAC = 1 - AUC). Screen MTS010 preferred when a drug-line pair has several (PRISM readme).
- Activity: VIPER (decoupler) with Bahl's full TCGA regulons on DepMap expression of the TNBC lines,
  genes z-scored across those lines.
- Statistic: Harrell's C between activity and sensitivity; two-sided permutation p (10,000 permutations).
- Replication: C on the reported side of 0.5 and p < 0.05. Primary outcome: replication rate with
  Clopper-Pearson 95% CI, vs 1,000 random sets of equal size (random PRISM drug x HC-TF pairs, each given
  the direction of a randomly chosen reported pair).
- Pre-specified: mTOR inhibitors x NFE2L3 and CBFB.
"""
import sys, re, numpy as np, pandas as pd, decoupler as dc
from scipy.stats import binomtest
root = sys.argv[1]; R = f'{root}/data/depmap_raw'
rng = np.random.default_rng(20261106)
norm = lambda s: re.sub(r'[^a-z0-9]', '', str(s).lower())

g = pd.read_csv(f'{root}/groups/cell_line_groups_FROZEN_2026-11-06.csv')
tnbc = set(g.loc[g.group == 'TNBC', 'ModelID'])
net = pd.read_csv(f'{root}/data/bahl2026/SuppData2_TCGA_VIPER_network_full.csv'); net.columns = ['source', 'target', 'weight']
net = net.groupby(['source', 'target'], as_index=False).weight.mean()
hc = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())

# Bahl pairs
s7 = pd.read_csv(f'{root}/data/bahl2026/SuppData7_drug_HCTF_concordance.csv')
s7 = s7[s7.FDR <= 0.05].copy(); s7['dkey'] = s7.drug.map(norm)
pairs = s7.groupby(['dkey', 'TR']).agg(drug=('drug', 'first'), n_cohorts=('Cohort', 'nunique'),
                                       mean_c=('c_index', 'mean'), cohorts=('Cohort', lambda x: '|'.join(sorted(set(x))))).reset_index()
pairs = pairs[pairs.n_cohorts >= 2].copy()
pairs['reported_dir'] = np.sign(pairs.mean_c - 0.5)
print(f'Bahl pairs significant in >= 2 datasets: {len(pairs)} (paper reports 176); TFs {pairs.TR.nunique()}, drugs {pairs.dkey.nunique()}')

# PRISM
p = pd.read_csv(f'{R}/prism20q2_secondary_dose_response_curve_parameters.csv',
                usecols=['depmap_id', 'screen_id', 'auc', 'name', 'moa', 'target', 'passed_str_profiling'])
p = p[(p.passed_str_profiling.astype(str).str.upper() == 'TRUE') & p.depmap_id.isin(tnbc) & p.auc.notna()].copy()
p['dkey'] = p.name.map(norm); p['pref'] = (p.screen_id == 'MTS010').astype(int)
p = p.sort_values('pref', ascending=False).groupby(['dkey', 'depmap_id'], as_index=False).agg(auc=('auc', 'first'), name=('name', 'first'), moa=('moa', 'first'))
sens = p.pivot(index='depmap_id', columns='dkey', values='auc').rsub(1)
names = p.groupby('dkey').name.first(); moas = p.groupby('dkey').moa.first()

# activity on TNBC lines
ex = pd.read_csv(f'{R}/depmap23q2_OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
ex.columns = [re.sub(r' \(\d+\)$', '', c) for c in ex.columns]; ex = ex.T.groupby(level=0).mean().T
lines = sorted(tnbc & set(ex.index) & set(sens.index))
X = ex.loc[lines]; X = X.loc[:, X.std() > 0]; Z = (X - X.mean()) / X.std()
res = dc.mt.viper(Z, net[net.target.isin(Z.columns)], tmin=5); nes = res[0] if isinstance(res, tuple) else res
print(f'TNBC lines with PRISM + expression: {len(lines)}')

I, J = None, None
def cidx_perm(x, y, B=10_000):
    ok = ~(np.isnan(x) | np.isnan(y)); x, y = x[ok], y[ok]; n = len(x)
    if n < 6: return np.nan, np.nan, n
    i, j = np.triu_indices(n, 1)
    def c(yy):
        dy = np.sign(yy[..., i] - yy[..., j]); dx = np.sign(x[i] - x[j]); s = dx * dy
        return ((s > 0) + 0.5 * ((dx == 0) & (dy != 0))).sum(-1) / (dy != 0).sum(-1)
    obs = c(y)
    perm = np.array([rng.permutation(y) for _ in range(B)])
    null = c(perm)
    pval = (1 + np.sum(np.abs(null - 0.5) >= abs(obs - 0.5) - 1e-12)) / (B + 1)
    return obs, pval, n

pairs['in_prism'] = pairs.dkey.isin(sens.columns) & pairs.TR.isin(nes.columns)
rows = []
for _, r in pairs[pairs.in_prism].iterrows():
    c, pv, n = cidx_perm(nes.loc[lines, r.TR].to_numpy(float), sens.loc[lines, r.dkey].to_numpy(float))
    rows.append(dict(drug=r.drug, prism_name=names[r.dkey], TF=r.TR, cohorts=r.cohorts, reported_mean_C=r.mean_c,
                     reported_dir=int(r.reported_dir), prism_C=c, prism_p=pv, n_lines=n,
                     replicated=bool(pv < 0.05 and np.sign(c - 0.5) == r.reported_dir)))
t = pd.DataFrame(rows); t.to_csv(f'{root}/results/10_h4_prism_pairs.csv', index=False, float_format='%.4g')
k, m = int(t.replicated.sum()), len(t)
ci = binomtest(k, m).proportion_ci(method='exact')
same_dir = int((np.sign(t.prism_C - 0.5) == t.reported_dir).sum())
print(f'Pairs testable in PRISM: {m} of {len(pairs)}; replicated: {k} ({k/m:.3f}; 95% CI {ci.low:.3f}-{ci.high:.3f}); same direction regardless of p: {same_dir}/{m}')

# null: random PRISM drug x HC-TF pairs
hcn = [x for x in hc if x in nes.columns]; dks = list(sens.columns); dirs = t.reported_dir.to_numpy()
null_rates = []
for b in range(1000):
    hits = 0
    for _ in range(m):
        tf, dk, d = rng.choice(hcn), rng.choice(dks), rng.choice(dirs)
        c, pv, n = cidx_perm(nes.loc[lines, tf].to_numpy(float), sens.loc[lines, dk].to_numpy(float), B=200)
        hits += int(pv < 0.05 and np.sign(c - 0.5) == d)
    null_rates.append(hits / m)
null_rates = np.array(null_rates)
p_emp = (1 + np.sum(null_rates >= k / m)) / 1001
print(f'Random-pair null replication rate: mean {null_rates.mean():.3f} (95% {np.percentile(null_rates,2.5):.3f}-{np.percentile(null_rates,97.5):.3f}); one-sided p (observed >= null) = {p_emp:.3f}')
pd.DataFrame([dict(n_pairs_reported=len(pairs), n_testable=m, n_replicated=k, rate=k/m, ci_lo=ci.low, ci_hi=ci.high,
                   same_direction=same_dir, null_mean=null_rates.mean(), null_lo=np.percentile(null_rates, 2.5),
                   null_hi=np.percentile(null_rates, 97.5), p_vs_null=p_emp, n_lines=len(lines))]).to_csv(f'{root}/results/10_h4_summary.csv', index=False)

# pre-specified mTOR x NFE2L3/CBFB
mt = [d for d in dks if isinstance(moas.get(d), str) and 'mtor' in moas[d].lower()]
rows = []
for d in mt:
    for tf in ['NFE2L3', 'CBFB']:
        if tf in nes.columns:
            c, pv, n = cidx_perm(nes.loc[lines, tf].to_numpy(float), sens.loc[lines, d].to_numpy(float))
            rows.append(dict(drug=names[d], moa=moas[d], TF=tf, C=c, p=pv, n=n))
mtor = pd.DataFrame(rows).sort_values(['TF', 'p']); mtor.to_csv(f'{root}/results/10_h4_mtor_prespecified.csv', index=False, float_format='%.4g')
pd.set_option('display.width', 200)
print('\nPre-specified mTOR inhibitors x NFE2L3 / CBFB (C > 0.5 = higher activity, more sensitive):')
print(mtor.to_string(index=False, float_format='%.3g'))
print('\nReplicated pairs:'); print(t[t.replicated][['drug', 'TF', 'cohorts', 'reported_mean_C', 'prism_C', 'prism_p', 'n_lines']].to_string(index=False, float_format='%.3g'))
