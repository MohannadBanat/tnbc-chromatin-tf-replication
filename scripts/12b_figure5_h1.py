"""Figure 5: cross-laboratory chromatin (H1) - Spearman correlations among the four shared lines in both labs."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
root = sys.argv[1]; F = f'{root}/figures'
INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'savefig.dpi': 300})
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('blue', ['#f0efec', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.7))
for ax, (fn, title, letter) in zip(axs, [('raw', 'All consensus peaks', 'A'),
                                         ('centred', 'Line-specific (centred within lab)', 'B')]):
    M = pd.read_csv(f'{root}/results/14_h1_shared_lines_corr_{fn}.csv', index_col=0)
    lab = [c.replace('_', ' ').replace('MDAMB', 'MDA-MB-').replace('BT549', 'BT-549') for c in M.columns]
    M.index = M.columns
    v = M.to_numpy(); vmin = np.floor(v[~np.eye(len(v), dtype=bool)].min() * 10) / 10
    vd = v.copy(); np.fill_diagonal(vd, np.nan); cm = cmap.copy(); cm.set_bad('white')
    im = ax.imshow(vd, cmap=cm, vmin=vmin, vmax=np.nanmax(vd))
    ax.set_xticks(range(len(lab))); ax.set_xticklabels(lab, rotation=55, ha='right', fontsize=7)
    ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab, fontsize=7)
    for i in range(len(v)):
        for j in range(len(v)):
            if i != j: ax.text(j, i, f'{v[i, j]:.2f}'.replace('-', '\u2212'), ha='center', va='center', fontsize=6.3, color='white' if v[i, j] > vmin + 0.65 * (np.nanmax(vd) - vmin) else INK)
    n = len(v) // 2
    ax.axhline(n - 0.5, color='white', lw=2); ax.axvline(n - 0.5, color='white', lw=2)
    ax.set_title(title); ax.text(-0.32, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight='bold')
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03); cb.set_label('Spearman ρ', fontsize=7); cb.ax.tick_params(labelsize=7)
    for s in ax.spines.values(): s.set_visible(False)
fig.tight_layout(); fig.savefig(f'{F}/Figure5_crosslab_chromatin.png'); print('ok')
