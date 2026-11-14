"""H1 (plan 2026-11-05): cross-laboratory chromatin reproducibility.
Bahl et al. ATAC-seq peaks (GSE323120; narrowPeak only, hg38) vs Lund ATAC-seq coverage (GSE337500; RPKM bigWig only, hg38).
The plan's fallback applies ("from the deposited bigWig signal at a common peak set"), because neither lab deposited reads/BAMs
in these series.
- Consensus peaks: Bahl peaks on chr1-22, X from all samples, merged; kept if called in >= 2 samples.
- Bahl value per region and sample: max signalValue (col 5) of its peaks inside the region, 0 if not called; mean over samples of a line.
- Lund value: mean bigWig coverage over the region per replicate (bases without data = 0), mean over replicates.
- Second region set: HC-TF promoters, TSS +/- 1 kb (NCBI RefSeq GRCh38 GTF, most upstream transcript start of the gene).
- Primary criterion (fixed in the plan): for each shared line X, Spearman(Lund X, Bahl X) > every within-lab different-line
  correlation involving X (Bahl X vs Bahl Y; Lund X vs Lund Y). Success if this holds for >= 3 of 4 shared lines.
- Exploratory: cross-lab identity - Spearman(Lund X, Bahl Y) for all Bahl lines Y; rank of the same line. Done on raw values and on
  values centred per region within each lab (line-specific accessibility), where the common open-chromatin component is removed.
- Not possible from deposited files (deviation): Jaccard of IDR peaks (Lund deposited no peaks), chromVAR (needs fragment counts).
"""
import sys, re, glob, gzip, os, numpy as np, pandas as pd, pyBigWig
from scipy.stats import spearmanr, rankdata
root = sys.argv[1]; G = f'{root}/data/geo'
CHR = [f'chr{i}' for i in range(1, 23)] + ['chrX']

def region_means(bw, chrs, starts, ends, span=10_000_000):
    """Mean coverage over [s, e) per region, bases without data counted as 0. Regions are grouped into windows of
    <= span bp; per window: dense per-base values and a float64 cumulative sum (bounded memory)."""
    out = np.zeros(len(chrs)); chrs = np.asarray(chrs); starts = np.asarray(starts, np.int64); ends = np.asarray(ends, np.int64)
    for c in np.unique(chrs):
        if c not in bw.chroms(): continue
        idx = np.where(chrs == c)[0]; idx = idx[np.argsort(starts[idx])]
        i = 0
        while i < len(idx):
            a = starts[idx[i]]; j = i
            while j < len(idx) and ends[idx[j]] - a <= span: j += 1
            j = max(j, i + 1); grp = idx[i:j]; b = ends[grp].max()
            v = bw.values(c, int(a), int(b), numpy=True); np.nan_to_num(v, copy=False)
            cs = np.concatenate([[0.0], np.cumsum(v, dtype=np.float64)])
            out[grp] = (cs[ends[grp] - a] - cs[starts[grp] - a]) / (ends[grp] - starts[grp]); i = j
    return out

import time; T0 = time.time()
def log(msg): print(f'[{time.time() - T0:6.0f}s] {msg}', flush=True)
def norm(s):
    s = re.sub(r'[^A-Z0-9]', '', s.upper())
    return {'MDA468': 'MDAMB468', 'HS578': 'HS578T', 'SUM149': 'SUM149PT', 'MB157': 'MDAMB157'}.get(s, s)
# ---- Bahl peaks
rows = []
for f in sorted(glob.glob(f'{G}/GSE323120/GSM*_peaks.filtered.merged.narrowPeak.gz')):
    b = os.path.basename(f); gsm = b.split('_')[0]
    line = norm(re.match(r'GSM\d+_(.+?)_\d+_[A-Z0-9]+_peaks', b).group(1))
    d = pd.read_csv(f, sep='\t', header=None, usecols=[0, 1, 2, 4], names=['chr', 's', 'e', 'sig'])
    d = d[d.chr.isin(CHR)]; d['sample'] = gsm; d['line'] = line; rows.append(d)
pk = pd.concat(rows, ignore_index=True)
smp = pk[['sample', 'line']].drop_duplicates()
print(f'Bahl: {len(smp)} samples, {smp.line.nunique()} lines, {len(pk)} peaks')
pk = pk.sort_values(['chr', 's']).reset_index(drop=True)
reg = np.zeros(len(pk), int); rid = -1
for c, sub in pk.groupby('chr', sort=False):
    idx = sub.index.to_numpy(); s = sub.s.to_numpy(); e = sub.e.to_numpy()
    newblk = np.r_[True, s[1:] > np.maximum.accumulate(e)[:-1]]
    reg[idx] = rid + np.cumsum(newblk); rid = reg[idx].max()
pk['reg'] = reg
R = pk.groupby('reg').agg(chr=('chr', 'first'), s=('s', 'min'), e=('e', 'max'), n=('sample', 'nunique'))
keep = R.index[R.n >= 2]; R = R.loc[keep]
print(f'Consensus regions (called in >= 2 samples): {len(R)} of {rid + 1} merged regions; median width {int((R.e - R.s).median())} bp')
S = pk[pk.reg.isin(keep)].groupby(['reg', 'sample']).sig.max().unstack(fill_value=0).reindex(keep, fill_value=0)
BL = S.T.groupby(smp.set_index('sample').line).mean().T                    # regions x Bahl lines
log('Bahl matrix built')
# ---- Lund coverage
lund = {}
for f in sorted(glob.glob(f'{G}/GSE337500/GSM*.bw')):
    m = re.match(r'GSM\d+_(.+)_rep(\d)\.bw', os.path.basename(f)); line = norm(m.group(1))
    bw = pyBigWig.open(f); vals = region_means(bw, R.chr.to_numpy(), R.s.to_numpy(), R.e.to_numpy())
    lund.setdefault(line, []).append(vals); bw.close(); log(f'Lund {os.path.basename(f)} done')
LU = pd.DataFrame({k: np.mean(v, 0) for k, v in lund.items()}, index=R.index)
rep_lund = {k: spearmanr(*v).statistic for k, v in lund.items() if len(v) == 2}
rep_bahl = {}
for l, sub in smp.groupby('line'):
    if len(sub) >= 2: rep_bahl[l] = spearmanr(S[sub['sample'].iloc[0]], S[sub['sample'].iloc[1]]).statistic
shared = sorted(set(BL.columns) & set(LU.columns))
print(f'Lund lines: {sorted(LU.columns)}; shared lines: {shared}')
print('Replicate agreement (Spearman over consensus regions): Lund', {k: round(v, 3) for k, v in rep_lund.items() if k in shared},
      '| Bahl', {k: round(v, 3) for k, v in rep_bahl.items() if k in shared})

def evaluate(B, L, label):
    rb = pd.DataFrame(rankdata(B, axis=0), index=B.index, columns=B.columns).corr()   # Spearman via ranks
    rl = pd.DataFrame(rankdata(L, axis=0), index=L.index, columns=L.columns).corr()
    out = []
    for x in shared:
        same = spearmanr(L[x], B[x]).statistic
        wb = rb.loc[x].drop(x); wl = rl.loc[x].drop(x)
        cross = pd.Series({y: spearmanr(L[x], B[y]).statistic for y in B.columns})
        out.append(dict(set=label, line=x, same_line_crosslab=same, max_within_Bahl=wb.max(), max_within_Bahl_line=wb.idxmax(),
                        max_within_Lund=wl.max(), max_within_Lund_line=wl.idxmax(), primary_pass=bool(same > max(wb.max(), wl.max())),
                        crosslab_rank_of_same_line=int((cross > same).sum() + 1), n_bahl_lines=len(cross),
                        best_crosslab_other=cross.drop(x).max(), best_crosslab_other_line=cross.drop(x).idxmax()))
    return pd.DataFrame(out)

res = [evaluate(BL, LU, 'consensus peaks, raw')]
lB = np.log1p(BL); lL = np.log1p(LU)
cB = lB.sub(lB.mean(1), axis=0); cL = lL.sub(lL.mean(1), axis=0)
res.append(evaluate(cB, cL, 'consensus peaks, centred within lab (exploratory)'))
# ---- HC-TF promoters
hc = sorted(pd.read_csv(f'{root}/data/bahl2026/SuppData5_TCGA_94HCTF_VIPER_network.csv')['Source (TR)'].unique())
acc = {f'NC_0000{i:02d}': f'chr{i}' for i in range(1, 23)}; acc['NC_000023'] = 'chrX'
tss = {}
with gzip.open(f'{G}/refseq/GRCh38_latest_genomic.gtf.gz', 'rt') as fh:
    for ln in fh:
        if ln.startswith('#'): continue
        p = ln.split('\t', 9)
        if p[2] != 'transcript': continue
        c = acc.get(p[0].split('.')[0])
        if c is None: continue
        g = re.search(r'gene_id "([^"]+)"', p[8]).group(1)
        if g not in hc: continue
        t = int(p[3]) if p[6] == '+' else int(p[4])
        if g not in tss: tss[g] = (c, t, p[6])
        else:
            c0, t0, st = tss[g]; tss[g] = (c0, min(t0, t) if st == '+' else max(t0, t), st)
P = pd.DataFrame([(g, c, t - 1000, t + 1000) for g, (c, t, _) in tss.items()], columns=['TF', 'chr', 's', 'e']).set_index('TF')
print(f'HC-TF promoters located: {len(P)} of {len(hc)}')
pkk = pk.set_index('chr')
pb = pd.DataFrame(0.0, index=P.index, columns=sorted(smp['sample']))
for tf, r in P.iterrows():
    sub = pk[(pk.chr == r.chr) & (pk.e > r.s) & (pk.s < r.e)]
    for smpl, v in sub.groupby('sample').sig.max().items(): pb.loc[tf, smpl] = v
PB = pb.T.groupby(smp.set_index('sample').line).mean().T
pl = {}
for f in sorted(glob.glob(f'{G}/GSE337500/GSM*.bw')):
    m = re.match(r'GSM\d+_(.+)_rep(\d)\.bw', os.path.basename(f)); line = norm(m.group(1)); bw = pyBigWig.open(f)
    pl.setdefault(line, []).append([float(np.nan_to_num(bw.values(r.chr, int(r.s), int(r.e), numpy=True)).mean()) for _, r in P.iterrows()]); bw.close()
PL = pd.DataFrame({k: np.mean(v, 0) for k, v in pl.items()}, index=P.index)
res.append(evaluate(PB, PL, 'HC-TF promoters (TSS +/- 1 kb)'))
out = pd.concat(res, ignore_index=True)
out.to_csv(f'{root}/results/14_h1_crosslab.csv', index=False, float_format='%.4g')
pd.set_option('display.width', 250); print(out.to_string(index=False, float_format='%.3f'))
for lab in out.set.unique():
    o = out[out.set == lab]; print(f'{lab}: primary criterion met for {o.primary_pass.sum()}/{len(o)} shared lines; same line ranked first cross-lab for {(o.crosslab_rank_of_same_line == 1).sum()}/{len(o)}')
pd.DataFrame(rankdata(np.c_[BL[shared], LU[shared]], axis=0), columns=[f'Bahl_{x}' for x in shared] + [f'Lund_{x}' for x in shared]).corr().to_csv(f'{root}/results/14_h1_shared_lines_corr_raw.csv', float_format='%.4g')
pd.DataFrame(np.c_[cB[shared], cL[shared]], columns=[f'Bahl_{x}' for x in shared] + [f'Lund_{x}' for x in shared]).corr('spearman').to_csv(f'{root}/results/14_h1_shared_lines_corr_centred.csv', float_format='%.4g')
pd.DataFrame([dict(line=k, lab='Lund', rep_rho=v) for k, v in rep_lund.items()] + [dict(line=k, lab='Bahl', rep_rho=v) for k, v in rep_bahl.items()]).to_csv(f'{root}/results/14_h1_replicates.csv', index=False, float_format='%.4g')
