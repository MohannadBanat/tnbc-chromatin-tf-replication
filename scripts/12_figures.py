"""Manuscript figures (300 dpi PNG). Colours: blue = this study, orange = reported by Bahl et al., grey = null/chance."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root = sys.argv[1]; F = f'{root}/figures'; import os; os.makedirs(F, exist_ok=True)
BLUE, ORANGE, GREY, INK, INK2 = '#2a78d6', '#eb6834', '#9a9893', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.titlesize': 9, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'savefig.dpi': 300})
def panel(ax, letter): ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom')

# ---------- Figure 1: PDX re-test ----------
r = pd.read_csv(f'{root}/results/01_everolimus_pdx_recomputed.csv')
s = r[r.FDR_reported <= 0.05].sort_values('p_reported', ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(5.4, 4.0))
y = np.arange(len(s))
for i, row in s.iterrows():
    ax.plot([row.p_reported, row.p_perm_at_reported_c], [i, i], color='#d6d5d0', lw=1.2, zorder=1)
ax.scatter(s.p_reported, y, s=28, color=ORANGE, zorder=3, label='Reported by Bahl et al.')
ax.scatter(s.p_perm_at_reported_c, y, s=28, color=BLUE, zorder=3, label='Permutation p at the reported C (this study)')
ax.set_xscale('log'); ax.set_yticks(y); ax.set_yticklabels(s.TF)
minp = 2 / 39916800
ax.axvline(minp, color=INK2, ls='--', lw=0.8); ax.text(minp * 1.4, len(s) - 0.4, 'smallest p any exact test\ncan give with n = 11', fontsize=7, color=INK2, va='top')
ax.axvline(0.05, color=GREY, ls=':', lw=0.8); ax.text(0.042, -0.95, 'p = 0.05', fontsize=7, color=INK2, ha='right', va='center'); ax.set_ylim(-1.35, len(s) - 0.45)
ax.set_xlabel('p-value (log scale)'); ax.set_xlim(1e-17, 2)
ax.set_title('Everolimus response, 11 TNBC PDXs: TFs reported significant')
ax.legend(frameon=False, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=2)
fig.tight_layout(); fig.savefig(f'{F}/Figure1_PDX_retest.png'); plt.close(fig)

# ---------- Figure 2: dependency (H3) ----------
e1 = pd.read_csv(f'{root}/results/05_h3_primary_h3a_enrichment.csv'); e2 = pd.read_csv(f'{root}/results/07ps_h3_primary_h3a_enrichment.csv')
sel = pd.read_csv(f'{root}/results/05_h3a_selectivity_per_TF.csv')
b = pd.read_csv(f'{root}/results/06_h3b_activity_dependency_per_TF.csv')
fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.0), gridspec_kw={'width_ratios': [1.0, 1.1, 1]})
ax = axs[0]; panel(ax, 'A')
labs = ['DepMap\nall', 'DepMap\nno CE', 'Score\nall', 'Score\nno CE']
rows = [e1.iloc[0], e1.iloc[1], e2.iloc[0], e2.iloc[1]]
for i, rr in enumerate(rows):
    ax.errorbar(i, rr.null_mean * 100, yerr=[[100 * (rr.null_mean - rr.null_lo)], [100 * (rr.null_hi - rr.null_mean)]], fmt='o', color=GREY, ms=5, capsize=3, lw=1)
    ax.scatter(i + 0.18, rr.rate * 100, s=34, color=BLUE, zorder=3)
ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=7); ax.set_xlim(-0.5, 3.6)
ax.set_ylabel('% of TFs that are TNBC dependencies'); ax.set_title('Enrichment'); ax.set_ylim(-0.8, 17)
ax.scatter([], [], color=BLUE, s=20, label='93 HC-TFs tested'); ax.errorbar([], [], yerr=[], fmt='o', color=GREY, label='Matched TFs (95% range)')
ax.legend(frameon=False, fontsize=7, loc='upper center')
ax = axs[1]; panel(ax, 'B')
ce = sel.common_essential
ax.scatter(sel.median_luminal[~ce], sel.median_TNBC[~ce], s=12, color=BLUE, alpha=0.8, label='HC-TF')
ax.scatter(sel.median_luminal[ce], sel.median_TNBC[ce], s=22, facecolor='white', edgecolor=INK, lw=1, label='Common essential')
lim = [-2.6, 0.6]; ax.plot(lim, lim, color='#d6d5d0', lw=0.8, zorder=0)
ax.axhline(-0.5, color=GREY, ls=':', lw=0.8); ax.axvline(-0.5, color=GREY, ls=':', lw=0.8)
off = {'NFE2L3': (-34, -10), 'SRF': (5, -3), 'HMGA1': (5, 1), 'ILF2': (5, -6), 'MYC': (5, -3), 'MED30': (5, -3)}
for _, rr in sel[ce].iterrows(): ax.annotate(rr.TF, (rr.median_luminal, rr.median_TNBC), xytext=off.get(rr.TF, (4, -2)), textcoords='offset points', fontsize=6.5, color=INK2)
for t, o in [('SOX9', (8, 6)), ('E2F3', (8, -9))]:
    rr = sel[sel.TF == t].iloc[0]; ax.annotate(t, (rr.median_luminal, rr.median_TNBC), xytext=o, textcoords='offset points', fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle='-', color=GREY, lw=0.5))
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel('Median gene effect, luminal lines'); ax.set_ylabel('Median gene effect, TNBC lines')
ax.set_title('TNBC vs luminal'); ax.legend(frameon=False, fontsize=7, loc='lower right')
ax = axs[2]; panel(ax, 'C')
bins = np.linspace(-0.8, 0.8, 33)
ax.hist(b[~b.HC].rho_activity.dropna(), bins=bins, density=True, color=GREY, alpha=0.55, label='Other TFs')
ax.hist(b[b.HC].rho_activity.dropna(), bins=bins, density=True, histtype='step', color=BLUE, lw=1.6, label='HC-TFs')
ax.axvline(0, color=INK2, lw=0.6)
ax.set_xlabel('Spearman ρ, activity vs effect'); ax.set_ylabel('Density'); ax.set_title('Activity vs dependency')
ax.set_ylim(0, ax.get_ylim()[1] * 1.25)
ax.legend(frameon=False, fontsize=7, loc='upper left')
fig.tight_layout(); fig.savefig(f'{F}/Figure2_dependency.png'); plt.close(fig)

# ---------- Figure 3: drug replication (H4) ----------
t = pd.read_csv(f'{root}/results/10_h4_prism_pairs.csv'); sm = pd.read_csv(f'{root}/results/10_h4_summary.csv').iloc[0]
# random-pair null from the full 10,000-permutation re-run (script 10b)
nb = pd.read_csv(f'{root}/results/10b_h4_null_10k_summary.csv').iloc[0]
for k in ['null_mean', 'null_lo', 'null_hi', 'p_vs_null']: sm[k] = nb[k]
fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.9), gridspec_kw={'width_ratios': [1.4, 1]})
ax = axs[0]; panel(ax, 'A')
rep = t.replicated
ax.scatter(t.reported_mean_C[~rep], t.prism_C[~rep], s=12, color=BLUE, alpha=0.7, label='Not replicated')
ax.scatter(t.reported_mean_C[rep], t.prism_C[rep], s=34, color=ORANGE, edgecolor='white', lw=0.8, label='Replicated (p < 0.05, same side)')
for _, rr in t[rep].iterrows(): ax.annotate(f'{rr.drug}–{rr.TF}', (rr.reported_mean_C, rr.prism_C), xytext=(4, 2), textcoords='offset points', fontsize=6.5, color=INK2)
ax.axhline(0.5, color=GREY, lw=0.6); ax.axvline(0.5, color=GREY, lw=0.6)
ax.set_xlabel('C reported by Bahl et al. (mean of 2–3 screens)'); ax.set_ylabel('C in PRISM (16 TNBC lines max)')
ax.set_title('148 drug–TF pairs re-tested in PRISM'); ax.legend(frameon=False, fontsize=7, loc='upper left')
ax = axs[1]; panel(ax, 'B')
ax.errorbar(0, sm.null_mean * 100, yerr=[[100 * (sm.null_mean - sm.null_lo)], [100 * (sm.null_hi - sm.null_mean)]], fmt='o', color=GREY, capsize=3, ms=5)
ax.errorbar(1, sm.rate * 100, yerr=[[100 * (sm.rate - sm.ci_lo)], [100 * (sm.ci_hi - sm.rate)]], fmt='o', color=BLUE, capsize=3, ms=6)
ax.set_xticks([0, 1]); ax.set_xticklabels(['Random drug–TF\npairs (95% range)', 'Bahl et al.\npairs (95% CI)']); ax.set_xlim(-0.6, 1.6)
ax.set_ylabel('% of pairs replicated'); ax.set_title('Replication rate (H4)')
ax.text(0.5, 0.92, f'p = {sm.p_vs_null:.2f}', transform=ax.transAxes, ha='center', fontsize=7, color=INK2)
fig.tight_layout(); fig.savefig(f'{F}/Figure3_drug_replication.png'); plt.close(fig)

# ---------- Figure 4: activity reproducibility (H2) ----------
h2 = pd.read_csv(f'{root}/results/11_h2_depmap_activity_agreement.csv', index_col=0)
fig, ax = plt.subplots(figsize=(3.4, 2.7))
bins = np.linspace(-0.2, 1, 25)
ax.hist(h2[~h2.HC].rho_released_vs_DepMap, bins=bins, color=GREY, alpha=0.55, label='Other TNBC-specific TFs')
ax.hist(h2[h2.HC].rho_released_vs_DepMap, bins=bins, histtype='step', color=BLUE, lw=1.6, label='HC-TFs')
ax.set_xlabel('Spearman ρ, released vs recomputed activity\n(49 shared cell lines)'); ax.set_ylabel('Number of TFs')
ax.set_title('TF activity is reproducible (H2)'); ax.legend(frameon=False, fontsize=7, loc='upper left')
fig.tight_layout(); fig.savefig(f'{F}/Figure4_activity_reproducibility.png'); plt.close(fig)
print('figures written:', sorted(os.listdir(F)))
