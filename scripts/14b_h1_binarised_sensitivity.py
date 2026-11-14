"""H1 sensitivity analysis (14b): binarised accessibility. Same regions and matrices as script 14.
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


# ---- cache matrices for reuse
import pickle
pickle.dump(dict(BL=BL, LU=LU, S=S, smp=smp, lund=lund), open(f'{root}/results/14b_h1_matrices.pkl', 'wb'))
log('matrices cached')
# ---- binarise. Bahl: called (> 0) per sample. Lund: per replicate, the top-k regions by coverage are 'called', with k = the
# median number of consensus regions called per Bahl sample, so both laboratories call the same number of regions.
called_B = (S > 0)
k = int(called_B.sum(0).median())
def topk(v):
    thr = np.sort(v)[-k]; return v >= thr
LB = called_B.T.groupby(smp.set_index('sample').line).mean().T >= 0.5           # line called if called in >= half its samples
LL = pd.DataFrame({ln: (np.mean([topk(v) for v in reps], 0) >= 0.5) for ln, reps in lund.items()}, index=R.index)
print(f'k (regions called per sample, Bahl median) = {k}; called per line: Bahl median {int(LB.sum().median())}, Lund median {int(LL.sum().median())}')
def jac(a, b): a = a.to_numpy(bool); b = b.to_numpy(bool); return (a & b).sum() / (a | b).sum()
rows = []
for x in shared:
    same = jac(LL[x], LB[x])
    wb = pd.Series({y: jac(LB[x], LB[y]) for y in LB.columns if y != x}); wl = pd.Series({y: jac(LL[x], LL[y]) for y in LL.columns if y != x})
    cross = pd.Series({y: jac(LL[x], LB[y]) for y in LB.columns})
    rows.append(dict(set='binarised, Jaccard', line=x, same_line_crosslab=same, max_within_Bahl=wb.max(), max_within_Bahl_line=wb.idxmax(),
                     max_within_Lund=wl.max(), max_within_Lund_line=wl.idxmax(), primary_pass=bool(same > max(wb.max(), wl.max())),
                     crosslab_rank_of_same_line=int((cross > same).sum() + 1), n_bahl_lines=len(cross),
                     best_crosslab_other=cross.drop(x).max(), best_crosslab_other_line=cross.drop(x).idxmax()))
out = pd.DataFrame(rows); out.to_csv(f'{root}/results/14b_h1_binarised.csv', index=False, float_format='%.4g')
pd.set_option('display.width', 250); print(out.to_string(index=False, float_format='%.3f'))
print(f'binarised: primary criterion met for {out.primary_pass.sum()}/{len(out)}; same line ranked first cross-lab for {(out.crosslab_rank_of_same_line == 1).sum()}/{len(out)}')
