"""H3c agreement between CRISPR screens and H3d RNAi (DEMETER2), amendment 01.
- Agreement: Spearman of per-TF median gene effect across TNBC lines screened in BOTH DepMap 23Q2 and
  Project Score (same lines), for HC-TFs and for all TFs.
- RNAi: HC-TF median DEMETER2 effect across frozen TNBC lines; count <= -0.5; agreement with CRISPR.
Descriptive (no hypothesis test beyond Spearman), as pre-specified.
"""
import sys, re, numpy as np, pandas as pd
from scipy.stats import spearmanr
root = sys.argv[1]; R = f'{root}/data/depmap_raw'
g = pd.read_csv(f'{root}/groups/cell_line_groups_FROZEN_2026-11-06.csv')
tnbc = set(g.loc[g.group == 'TNBC', 'ModelID'])
universe = set(pd.read_csv(f'{root}/data/bahl2026/SuppData2_TCGA_VIPER_network_full.csv').iloc[:, 0])
hc = set(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'])
clean = lambda c: re.sub(r' \([\d&]+\)$', '', c)
def load(f):
    d = pd.read_csv(f'{R}/{f}', index_col=0); d.columns = [clean(c) for c in d.columns]
    return d.loc[:, [c for c in d.columns if c in universe]].T.groupby(level=0).mean().T
dm, ps = load('depmap23q2_CRISPRGeneEffect.csv'), load('projectscore_chronos_gene_effect.csv')
both = sorted(tnbc & set(dm.index) & set(ps.index))
genes = sorted(set(dm.columns) & set(ps.columns))
a, b = dm.loc[both, genes].median(), ps.loc[both, genes].median()
rho_all = spearmanr(a, b).statistic; hcg = [x for x in genes if x in hc]; rho_hc = spearmanr(a[hcg], b[hcg]).statistic
# per-line agreement too: per TF across shared lines
per_tf = pd.Series({t: spearmanr(dm.loc[both, t], ps.loc[both, t]).statistic for t in hcg})
print(f'TNBC lines in both CRISPR screens: {len(both)} ({", ".join(both)})')
print(f'Agreement of per-TF median effect: all TFs rho = {rho_all:.3f} (n = {len(genes)}); HC-TFs rho = {rho_hc:.3f} (n = {len(hcg)})')
print(f'Per-HC-TF across-line agreement: median rho = {per_tf.median():.3f}')

d2 = pd.read_csv(f'{R}/demeter2_gene_effect.csv', index_col=0)
if not str(d2.index[0]).startswith('ACH') and not str(d2.index[0]).endswith(tuple('ABCDEFGHIJKLMNOPQRSTUVWXYZ')):
    pass
info = pd.read_csv(f'{R}/demeter2_sample_info.csv', index_col=0)
first = str(d2.index[0])
if '(' in first:            # genes in rows
    d2 = d2.T
d2.columns = [clean(c) for c in d2.columns]
if not str(d2.index[0]).startswith('ACH'):
    d2.index = d2.index.map(dict(zip(info.CCLE_ID, info.index)))
d2 = d2.loc[:, [c for c in d2.columns if c in universe]].T.groupby(level=0).mean().T
tl = sorted(tnbc & set(d2.index.dropna()))
hc2 = [x for x in d2.columns if x in hc]
m2 = d2.loc[tl, hc2].median()
cr = dm.loc[[x for x in tl if x in dm.index], [x for x in hc2 if x in dm.columns]].median() if any(x in dm.index for x in tl) else None
print(f'\nRNAi: TNBC lines {len(tl)}; HC-TFs measured {len(hc2)}; HC-TFs with median <= -0.5: {sorted(m2[m2 <= -0.5].index)}')
if cr is not None:
    common = [x for x in hc2 if x in cr.index]
    print(f'CRISPR vs RNAi per-HC-TF median effect (TNBC lines in both): rho = {spearmanr(cr[common], m2[common]).statistic:.3f} (n = {len(common)})')
out = pd.DataFrame({'median_DepMap': a[hcg], 'median_ProjectScore': b[hcg], 'per_TF_rho_lines': per_tf}).join(m2.rename('median_RNAi'), how='outer')
out.to_csv(f'{root}/results/08_h3c_h3d_agreement.csv', float_format='%.4g')
