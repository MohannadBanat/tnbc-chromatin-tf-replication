"""H4, second planned screen (plan 2026-11-05 section 2, H4: "PRISM and/or GDSC1"): replication of Bahl et al.'s
cross-dataset drug-TF pairs in GDSC1.

Added 2026-11-14 after co-author review (see DEVIATIONS.md). All rules below were fixed before any GDSC1
response value was related to TF activity; they mirror script 10 / 10b (PRISM).

Data: DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv.gz (Sanger fitted dose-response, release of 2020-06-02;
      figshare file 34355645 in doi:10.6084/m9.figshare.19345397). Rows with dataset == 'GDSC1' only.
- Lines: frozen TNBC group (groups/cell_line_groups_FROZEN_2026-11-06.csv), matched on BROAD_ID = DepMap ModelID,
  with DepMap 23Q2 expression.
- Drug matching: normalised drug_name, or any normalised entry of drug_synonyms, equals Bahl's normalised drug name.
  If a drug has several GDSC1 drug_ids, the drug_id with most TNBC lines is used (ties: lowest drug_id).
- Sensitivity = 1 - AUC (Bahl et al. use AAC = 1 - AUC).
- Activity: VIPER (decoupler) with Bahl's full TCGA regulons on DepMap expression of the TNBC lines with GDSC1 data,
  genes z-scored across those lines.
- Statistic: Harrell's C between activity and sensitivity; two-sided Monte Carlo permutation p (10,000 permutations,
  cell line as unit, p = (b+1)/(B+1)).
- Replication: C on the reported side of 0.5 and p < 0.05. Primary outcome: replication rate with Clopper-Pearson
  95% CI vs 1,000 random sets of equal size (random GDSC1 drug x HC-TF pairs, direction drawn from the reported
  pairs), 10,000 permutations per null pair.
- Pre-specified: mTOR inhibitors x NFE2L3 and CBFB (putative_target contains MTOR).
- Descriptive split (not pre-specified): pairs whose reported significance included Bahl's GDSC1000 cohort versus
  pairs significant in CTRPv2 + GRAY only, because GDSC1 is the same Sanger screening platform as GDSC1000.
"""
import sys, re, numpy as np, pandas as pd, decoupler as dc
from scipy.stats import binomtest
root = sys.argv[1]; R = f'{root}/data/depmap_raw'
rng = np.random.default_rng(20261114)
norm = lambda s: re.sub(r'[^a-z0-9]', '', str(s).lower())

g = pd.read_csv(f'{root}/groups/cell_line_groups_FROZEN_2026-11-06.csv')
tnbc = set(g.loc[g.group == 'TNBC', 'ModelID'])
net = pd.read_csv(f'{root}/data/bahl2026/SuppData2_TCGA_VIPER_network_full.csv'); net.columns = ['source', 'target', 'weight']
net = net.groupby(['source', 'target'], as_index=False).weight.mean()
hc = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())

s7 = pd.read_csv(f'{root}/data/bahl2026/SuppData7_drug_HCTF_concordance.csv')
s7 = s7[s7.FDR <= 0.05].copy(); s7['dkey'] = s7.drug.map(norm)
pairs = s7.groupby(['dkey', 'TR']).agg(drug=('drug', 'first'), n_cohorts=('Cohort', 'nunique'), mean_c=('c_index', 'mean'),
                                       cohorts=('Cohort', lambda x: '|'.join(sorted(set(x))))).reset_index()
pairs = pairs[pairs.n_cohorts >= 2].copy(); pairs['reported_dir'] = np.sign(pairs.mean_c - 0.5)
print(f'Bahl pairs significant in >= 2 datasets: {len(pairs)}')

d = pd.read_csv(f'{root}/data/gdsc_raw/DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv.gz', low_memory=False,
                usecols=['dataset', 'drug_id', 'drug_name', 'drug_synonyms', 'putative_target', 'BROAD_ID', 'AUC'])
d = d[(d.dataset == 'GDSC1') & d.BROAD_ID.isin(tnbc) & d.AUC.notna()].copy()
print(f'GDSC1 rows in frozen TNBC lines: {len(d)}; lines {d.BROAD_ID.nunique()}; drug_ids {d.drug_id.nunique()}')
d = d.groupby(['drug_id', 'BROAD_ID'], as_index=False).agg(AUC=('AUC', 'mean'), drug_name=('drug_name', 'first'),
                                                           drug_synonyms=('drug_synonyms', 'first'), putative_target=('putative_target', 'first'))
meta = d.groupby('drug_id').agg(drug_name=('drug_name', 'first'), syn=('drug_synonyms', 'first'),
                                target=('putative_target', 'first'), n=('BROAD_ID', 'nunique')).reset_index()
sens_id = d.pivot(index='BROAD_ID', columns='drug_id', values='AUC').rsub(1)
key2id = {}
for _, r in meta.sort_values(['n', 'drug_id'], ascending=[False, True]).iterrows():
    keys = {norm(r.drug_name)} | ({norm(x) for x in str(r.syn).split(',')} if isinstance(r.syn, str) else set())
    for k in keys - {''}: key2id.setdefault(k, r.drug_id)
name_of = meta.set_index('drug_id').drug_name; target_of = meta.set_index('drug_id').target

ex = pd.read_csv(f'{R}/depmap23q2_OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
ex.columns = [re.sub(r' \(\d+\)$', '', c) for c in ex.columns]; ex = ex.T.groupby(level=0).mean().T
lines = sorted(tnbc & set(ex.index) & set(sens_id.index))
X = ex.loc[lines]; X = X.loc[:, X.std() > 0]; Z = (X - X.mean()) / X.std()
res = dc.mt.viper(Z, net[net.target.isin(Z.columns)], tmin=5); nes = res[0] if isinstance(res, tuple) else res
print(f'TNBC lines with GDSC1 + expression: {len(lines)}')

def cidx_perm(x, y, gen, B=10_000):
    ok = ~(np.isnan(x) | np.isnan(y)); x, y = x[ok], y[ok]; n = len(x)
    if n < 6: return np.nan, np.nan, n
    i, j = np.triu_indices(n, 1); dx = np.sign(x[i] - x[j])
    def c(yy):
        dy = np.sign(yy[..., i] - yy[..., j]); s = dx * dy
        return ((s > 0) + 0.5 * ((dx == 0) & (dy != 0))).sum(-1) / (dy != 0).sum(-1)
    obs = c(y); null = c(gen.permuted(np.tile(y, (B, 1)), axis=1))
    return obs, (1 + np.sum(np.abs(null - 0.5) >= abs(obs - 0.5) - 1e-12)) / (B + 1), n

pairs['gid'] = pairs.dkey.map(key2id)
pairs['in_gdsc1'] = pairs.gid.notna() & pairs.TR.isin(nes.columns)
rows = []
for _, r in pairs[pairs.in_gdsc1].iterrows():
    c, pv, n = cidx_perm(nes.loc[lines, r.TR].to_numpy(float), sens_id.loc[lines, r.gid].to_numpy(float), rng)
    rows.append(dict(drug=r.drug, gdsc1_name=name_of[r.gid], gdsc1_drug_id=int(r.gid), TF=r.TR, cohorts=r.cohorts,
                     reported_in_GDSC1000='GDSC1000' in r.cohorts, reported_mean_C=r.mean_c, reported_dir=int(r.reported_dir),
                     gdsc1_C=c, gdsc1_p=pv, n_lines=n))
t = pd.DataFrame(rows); t = t[t.gdsc1_C.notna()].copy()
t['replicated'] = (t.gdsc1_p < 0.05) & (np.sign(t.gdsc1_C - 0.5) == t.reported_dir)
t.to_csv(f'{root}/results/16_h4_gdsc1_pairs.csv', index=False, float_format='%.4g')
def rate(tt):
    k, m = int(tt.replicated.sum()), len(tt)
    if m == 0: return k, m, np.nan, np.nan
    ci = binomtest(k, m).proportion_ci(method='exact'); return k, m, ci.low, ci.high
k, m, lo, hi = rate(t)
same_dir = int((np.sign(t.gdsc1_C - 0.5) == t.reported_dir).sum())
print(f'Pairs testable in GDSC1: {m} of {len(pairs)}; replicated {k} ({k/m:.3f}; 95% CI {lo:.3f}-{hi:.3f}); same direction {same_dir}/{m}; n lines median {t.n_lines.median():.0f} (range {t.n_lines.min()}-{t.n_lines.max()})')
split = {}
for lab, tt in [('reported_in_GDSC1000', t[t.reported_in_GDSC1000]), ('CTRPv2_GRAY_only', t[~t.reported_in_GDSC1000])]:
    split[lab] = rate(tt); print(lab, split[lab])

hcn = [x for x in hc if x in nes.columns]; ids = list(sens_id.columns); dirs = t.reported_dir.to_numpy()
# The null is run in two parallel parts of 500 sets (mode '1' and '2', own generators), then combined (mode 'combine').
mode = sys.argv[2] if len(sys.argv) > 2 else 'combine'
if mode in ('1', '2'):
    rng2 = np.random.default_rng(int(f'20261114{mode}')); null_rates = []
    for b in range(500):
        hits = 0; done = 0
        while done < m:
            tf, gid, dd = rng2.choice(hcn), rng2.choice(ids), rng2.choice(dirs)
            c, pv, n = cidx_perm(nes.loc[lines, tf].to_numpy(float), sens_id.loc[lines, gid].to_numpy(float), rng2)
            if np.isnan(c): continue
            done += 1; hits += int(pv < 0.05 and np.sign(c - 0.5) == dd)
        null_rates.append(hits / m)
        if b % 50 == 49: print(f'part {mode}: {b+1} null sets done', flush=True)
    np.save(f'{root}/results/16_h4_gdsc1_null_part{mode}.npy', np.array(null_rates)); sys.exit(0)
null_rates = np.concatenate([np.load(f'{root}/results/16_h4_gdsc1_null_part{i}.npy') for i in '12']); assert len(null_rates) == 1000
p_emp = (1 + np.sum(null_rates >= k / m)) / 1001
print(f'Random-pair null (1,000 sets; 10,000 permutations per pair): mean {null_rates.mean():.4f} (95% {np.percentile(null_rates,2.5):.4f}-{np.percentile(null_rates,97.5):.4f}); one-sided p = {p_emp:.3f}')
pd.DataFrame([dict(n_pairs_reported=len(pairs), n_testable=m, n_replicated=k, rate=k/m, ci_lo=lo, ci_hi=hi, same_direction=same_dir,
                   null_mean=null_rates.mean(), null_lo=np.percentile(null_rates, 2.5), null_hi=np.percentile(null_rates, 97.5),
                   p_vs_null=p_emp, n_lines=len(lines), n_lines_median_per_pair=t.n_lines.median(),
                   k_in_GDSC1000=split['reported_in_GDSC1000'][0], m_in_GDSC1000=split['reported_in_GDSC1000'][1],
                   k_CTRPv2_GRAY_only=split['CTRPv2_GRAY_only'][0], m_CTRPv2_GRAY_only=split['CTRPv2_GRAY_only'][1])]
             ).to_csv(f'{root}/results/16_h4_gdsc1_summary.csv', index=False)

mt = [i for i in ids if isinstance(target_of.get(i), str) and 'MTOR' in target_of[i].upper()]
rows = []
for i in mt:
    for tf in ['NFE2L3', 'CBFB']:
        c, pv, n = cidx_perm(nes.loc[lines, tf].to_numpy(float), sens_id.loc[lines, i].to_numpy(float), rng)
        rows.append(dict(drug=name_of[i], drug_id=int(i), target=target_of[i], TF=tf, C=c, p=pv, n=n))
mtor = pd.DataFrame(rows).dropna(subset=['C']).sort_values(['TF', 'p'])
from statsmodels.stats.multitest import multipletests
mtor['FDR'] = multipletests(mtor.p, method='fdr_bh')[1]
mtor.to_csv(f'{root}/results/16_h4_gdsc1_mtor_prespecified.csv', index=False, float_format='%.4g')
pd.set_option('display.width', 220)
print('\nmTOR-targeting compounds x NFE2L3 / CBFB:'); print(mtor.to_string(index=False, float_format='%.3g'))
print('\nReplicated pairs:'); print(t[t.replicated][['drug', 'TF', 'cohorts', 'reported_mean_C', 'gdsc1_C', 'gdsc1_p', 'n_lines']].to_string(index=False, float_format='%.3g'))
