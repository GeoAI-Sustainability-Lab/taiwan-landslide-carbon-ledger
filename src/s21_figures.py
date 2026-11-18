# -*- coding: utf-8 -*-
"""s21: Figures 1, 2, 3, 5, 6 and 8 and Figures S1-S4 (English labels). Figure 2 is a schematic; all other panels are computed from the outputs in v2/."""
import json, numpy as np, pandas as pd, shapefile
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import matplotlib.patheffects as pe, pyproj
TOWGS = pyproj.Transformer.from_crs(3826, 4326, always_xy=True)
D = 'v2'; FIG = 'figures'; import os; os.makedirs(FIG, exist_ok=True)
from config import SRC
INK = '#222222'; MUTE = '#666666'; GRID = '#E3E3E3'
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.6, 'axes.edgecolor': INK, 'axes.labelcolor': INK, 'text.color': INK,
                     'xtick.color': INK, 'ytick.color': INK, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'legend.fontsize': 6.8,
                     'axes.grid': False, 'savefig.dpi': 300})
COL = {'BL': '#009E73', 'CF': '#0072B2', 'MX': '#56B4E9', 'BB': '#E69F00', 'gain': '#CC79A7', 'official': '#D55E00'}
LAB = {'BL': 'Broadleaf', 'CF': 'Conifer', 'MX': 'Conifer-broadleaf mixed', 'BB': 'Bamboo and bamboo mixed'}; GS = ['BL', 'CF', 'MX', 'BB']; MK = {'BL': 'o', 'CF': 's', 'MX': '^', 'BB': 'D'}
S = json.load(open(D + '/ledger_summary.json')); P = json.load(open(D + '/parameters.json')); Q = json.load(open(D + '/support_stats.json'))
T = pd.read_csv(D + '/ledger_annual.csv'); years = T.year.values
def save(fig, name):
    fig.savefig(f'{FIG}/{name}.png', dpi=300, bbox_inches='tight', facecolor='white'); fig.savefig(f'{FIG}/{name}.pdf', bbox_inches='tight', facecolor='white'); plt.close(fig)
def panel(ax, s, x=-0.13, y=1.03): ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
def hgrid(ax, axis='y'): ax.grid(True, axis=axis, color=GRID, lw=0.5, ls='-'); ax.set_axisbelow(True)

# ---------------- Fig. 1: maps and annual catalogued area ----------------
sf = shapefile.Reader(SRC + '/Reference_Boundaries/NaturalEarth_admin_10m/ne_10m_admin_0_countries.shp', encoding='latin1')
fld = [f[0] for f in sf.fields[1:]]; ia = fld.index('ADM0_A3'); tw = [s for s, r in zip(sf.shapes(), sf.records()) if r[ia] == 'TWN'][0]
def outline(a_):
    pts = np.array(tw.points); parts = list(tw.parts) + [len(pts)]
    for k in range(len(parts) - 1): a_.plot(pts[parts[k]:parts[k + 1], 0], pts[parts[k]:parts[k + 1], 1], color=INK, lw=0.6, zorder=5)
mp = pd.read_csv(D + '/map_loss_0p05deg.csv')
fig = plt.figure(figsize=(7.2, 6.9)); gs = fig.add_gridspec(2, 2, height_ratios=[2.35, 1], hspace=0.22, wspace=0.12)
for j, (colm, cmap, lab, vmin, vmax, fmt) in enumerate([('ha', 'Greys', 'Catalogued landslide area 2004–2025 (ha per 0.05° cell)', 3, 3000, None), ('kt', 'YlOrBr', 'Living-biomass loss 2004–2025 (kt CO$_2$ per 0.05° cell)', 1, 1000, None)]):
    ax = fig.add_subplot(gs[0, j]); d = mp.assign(kt=mp.loss_t / 1e3); d = d[d[colm] >= vmin * 0.3]
    sc = ax.scatter(d.gx, d.gy, c=np.clip(d[colm], vmin, vmax), s=7.2, marker='s', cmap=cmap, norm=LogNorm(vmin, vmax), linewidths=0, zorder=3)
    outline(ax); ax.set_xlim(119.95, 122.05); ax.set_ylim(21.8, 25.4); ax.set_aspect(1 / np.cos(np.radians(23.6)))
    ax.set_xticks([120, 121, 122]); ax.set_yticks([22, 23, 24, 25]); ax.grid(True, color='#BBBBBB', lw=0.4, ls='-', zorder=1); ax.set_axisbelow(True)
    ax.set_xticklabels(['120°E', '121°E', '122°E']); ax.set_yticklabels(['22°N', '23°N', '24°N', '25°N'] if j == 0 else [])
    for sp in ax.spines.values(): sp.set_visible(True)
    cb = fig.colorbar(sc, ax=ax, orientation='horizontal', pad=0.06, fraction=0.04, aspect=32); cb.set_label(lab, fontsize=7); cb.outline.set_linewidth(0.5)
    panel(ax, 'ab'[j], x=-0.02 if j else -0.12, y=1.01)
    # 100 km scale bar
    x0s = 120.03; y0s = 25.18; ln = 100 / (111.32 * np.cos(np.radians(22))); ax.plot([x0s, x0s + ln], [y0s, y0s], color=INK, lw=1.2, zorder=6); ax.text(x0s + ln / 2, y0s + 0.05, '100 km', ha='center', fontsize=6.5, zorder=6)
    if j == 0:   # windows of the two case sites (Figs 4 and 7)
        for (r0_, r1_, c0_, c1_), lab_, dx_, ha_ in (((3650, 3900, 5380, 5640), 'Fig. 4', 0.09, 'left'), ((7500, 7800, 2300, 2700), 'Fig. 7', -0.11, 'right')):
            xs_ = 151050.0 + np.array([c0_, c1_, c1_, c0_, c0_]) * 30; ys_ = 2798910.0 - np.array([r0_, r0_, r1_, r1_, r0_]) * 30; lo_, la_ = TOWGS.transform(xs_, ys_)
            ax.plot(lo_, la_, color=INK, lw=0.9, zorder=7, path_effects=[pe.Stroke(linewidth=2.0, foreground='white'), pe.Normal()])
            ax.text(lo_.max() + dx_ if dx_ > 0 else lo_.min() + dx_, la_.mean(), lab_, ha=ha_, va='center', fontsize=6.5, zorder=7, path_effects=[pe.Stroke(linewidth=2.0, foreground='white'), pe.Normal()])
ax = fig.add_subplot(gs[1, :]); B = pd.read_csv(D + '/ledger_breakdown.csv'); tc = B[B.dim == 'treecover2000'].pivot_table(index='year', columns='class', values='area_ha', aggfunc='sum').reindex(years).fillna(0)
shades = ['#E0E0E0', '#BDBDBD', '#969696', '#636363', '#252525']; bot = np.zeros(len(years))
for c_, sh in zip(['0-25', '25-50', '50-75', '75-90', '90-100'], shades):
    ax.bar(years, tc[c_].values / 1e3, bottom=bot, color=sh, width=0.78, edgecolor='white', linewidth=0.5, label=c_.replace('-', '–') + '%'); bot += tc[c_].values / 1e3
for y in (2007, 2009, 2016, 2024): ax.text(y, bot[y - 2004] + 0.6, f'{bot[y - 2004]:.1f}', ha='center', fontsize=6.5)
ax.set_ylabel('Catalogued landslide area\n(10$^3$ ha yr$^{-1}$)'); ax.set_xticks(years[::2]); ax.set_xlim(2003.3, 2025.7); ax.set_ylim(0, 42); hgrid(ax)
ax.legend(title='Tree cover in 2000', ncol=5, loc='upper right', title_fontsize=6.8, columnspacing=1.0, handlelength=1.2); panel(ax, 'c', x=-0.07, y=1.0)
save(fig, 'Fig1_maps_area')

# ---------------- Fig. 2: workflow (schematic; nouns only; box size = role; colours = data origin) ----------------
fig, ax = plt.subplots(figsize=(7.1, 4.35)); ax.set_xlim(-0.5, 127.5); ax.set_ylim(6.5, 79.5); ax.axis('off')
C_IN, C_NIR, C_CHK, C_DER = '#E4E4E4', '#F7D3BA', '#ECECCF', '#FFFFFF'          # input data, national inventory (vermillion tint, as in Fig. 5), independent check data (olive tint), derived quantity
E_CHK, E_NIR = '#7A7A14', '#D55E00'                                              # olive (independent checks only) and vermillion (national inventory, as in Fig. 5)
def box(xc, yc, w, h, txt, fc, lw=0.7, fs=5.9, ec=INK, bold=False, ls='-'):
    ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle='round,pad=0.15,rounding_size=0.9', fc=fc, ec=ec, lw=lw, ls=ls, zorder=3))
    ax.text(xc, yc, txt, ha='center', va='center', fontsize=fs, linespacing=1.35, zorder=4, fontweight='bold' if bold else 'normal')
def arr(p0, p1, dashed=False):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=7.5, lw=0.75 if not dashed else 0.9, color=INK if not dashed else E_CHK,
                                 linestyle=(0, (3.2, 2.2)) if dashed else '-', shrinkA=0, shrinkB=0.6, zorder=2))
X = [12.5, 38, 63.5, 89, 115.5]; W, H = 23.5, 8.6; WC, HC = 22.0, 7.2; R1, R2, R3, R4, R5 = 73.5, 58.5, 43.5, 28.5, 13.5
top = lambda y, h=H: y + h / 2; bot = lambda y, h=H: y - h / 2; DSH = (0, (3.2, 1.6))
# row 1: input data
for x, t in zip(X[:4], ['Event-based landslide\ncatalogue 2004–2025', 'Forest inventory maps\nvolume and forest type', 'Tree cover 2000\nGlobal Forest Change', 'Canopy height CHMv2\nimage years 2010–2019']):
    box(x, R1, W, H, t, C_IN)
# row 2: one 30 m grid that carries every cell-level quantity
GX0, GX1, GH = X[0] - W / 2, X[3] + W / 2, 9.6
box((GX0 + GX1) / 2, R2, GX1 - GX0, GH, '30 m grid cells: landslide area A, failure date, stand stock S,\nforest type g, tree-cover class c, relative canopy height x', C_DER, fs=6.2)
for x in X[:4]: arr((x, bot(R1)), (x, R2 + GH / 2))
box(X[4], R2, WC, HC, 'GEDI L4A footprints\n2019–2025', C_CHK, ec=E_CHK, lw=0.9, ls=DSH); arr((X[4] - WC / 2, R2 - 2.2), (GX1, R2 - 2.2), dashed=True)
# row 3: parameters estimated from the grid
box(X[2], R3, W, H, 'Pre-failure biomass\nratio ρ', C_DER); box(X[3], R3, W, H, 'Recovery rate k and\nbare-scar baseline x$_0$', C_DER)
arr((X[2], R2 - GH / 2), (X[2], top(R3))); arr((X[3], R2 - GH / 2), (X[3], top(R3)))
# row 4: ledger (largest box) and the checks of rho and k
LX0, LX1 = X[0] - W / 2, X[1] + W / 2; LH = 12.0; LXC = (LX0 + LX1) / 2
box(LXC, R4, LX1 - LX0, LH, 'Annual living-biomass ledger 2004–2025\nloss L · re-accumulation G\nunrecovered stock D', C_DER, lw=1.5, fs=6.7)
arr((LXC, R2 - GH / 2), (LXC, R4 + LH / 2))
arr((X[2] - 7.5, bot(R3)), (LX1 - 3.5, R4 + LH / 2)); arr((X[3] - 7.5, bot(R3)), (LX1, R4 + 4.4))
box(X[2], R4, WC, HC, 'Third-inventory photo\npoints 1990–1994', C_CHK, ec=E_CHK, lw=0.9, ls=DSH); box(X[3], R4, WC, HC, 'Annual landslide\nmaps 2017–2024', C_CHK, ec=E_CHK, lw=0.9, ls=DSH)
arr((X[2], R4 + HC / 2), (X[2], bot(R3)), dashed=True); arr((X[3], R4 + HC / 2), (X[3], bot(R3)), dashed=True)
# row 5: outputs and the national inventory
box(X[0], R5, W, H, '95% intervals and\n12 alternative settings', C_DER); box(X[1], R5, W, H, 'National greenhouse\ngas inventory 2026', C_NIR, ec=E_NIR, lw=0.9); box(X[2], R5, W, H, 'Inventory comparison', C_DER)
arr((X[0], R4 - LH / 2), (X[0], top(R5))); arr((LX1 - 3.5, R4 - LH / 2), (X[2] - 7.5, top(R5)))
arr((X[1] - W / 2, R5), (X[0] + W / 2, R5)); arr((X[1] + W / 2, R5), (X[2] - W / 2, R5))
# legend (colour = origin of the box content; line = relation)
lg = [Patch(fc=C_IN, ec=INK, lw=0.7, label='Input data'), Patch(fc=C_NIR, ec=E_NIR, lw=0.9, label='National inventory data'), Patch(fc=C_CHK, ec=E_CHK, lw=0.9, ls=(0, (3.2, 1.6)), label='Independent data for checks'),
      Patch(fc=C_DER, ec=INK, lw=0.7, label='Derived quantity'), Patch(fc=C_DER, ec=INK, lw=1.5, label='Ledger, main output'),
      Line2D([], [], color=INK, lw=0.75, label='Calculation'), Line2D([], [], color=E_CHK, lw=0.9, ls=(0, (3.2, 2.2)), label='Independent check')]
ax.legend(handles=lg, loc='lower left', bbox_to_anchor=(0.835, 0.03), fontsize=6.0, handlelength=1.8, handleheight=1.1, labelspacing=0.75, borderaxespad=0)
save(fig, 'Fig2_workflow_schematic')

# ---------------- Fig. 3: pre-failure relative biomass ----------------
fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.5), sharey=True, gridspec_kw={'wspace': 0.08}); xs = np.arange(5); xl = ['0-25', '25-50', '50-75', '75-90', '>90']
for j, g in enumerate(GS):
    ax = axs[j]; r = np.array(P['rho'][g]); ci = np.array(P['rho_ci'][g]); rr = np.array(P['rho_ref'][g]); n = P['n_pre'][g]
    ax.plot(xs, rr, '-', color='#888888', lw=1.0, marker=MK[g], ms=4, mfc='white', mec='#888888', mew=0.8, label='Never-catalogued forest', zorder=2)
    ax.errorbar(xs, r, yerr=[r - ci[:, 0], ci[:, 1] - r], fmt='-', color=COL[g], lw=1.2, marker=MK[g], ms=4.5, mec='white', mew=0.6, capsize=2, elinewidth=0.8, label='Cells before failure', zorder=3)
    ax.axhline(1, color=INK, lw=0.5); ax.set_xticks(xs); ax.set_xticklabels(xl, fontsize=6.5); ax.set_title(LAB[g], fontsize=7.5); ax.set_ylim(0, 1.3); hgrid(ax)
    ax.text(0.03, 0.965, f"n = {sum(n):,}".replace(",", "\u2009"), transform=ax.transAxes, fontsize=6.3, va='top', color=MUTE)
    if j == 0: ax.set_ylabel('Relative biomass before failure, $\\rho$\n(local reference forest = 1)')
fig.supxlabel('Tree cover in 2000 (%)', fontsize=7.5, y=-0.06)
fig.legend(handles=[Line2D([], [], color=INK, marker='o', ms=4.5, mec='white', mew=0.6, lw=1.2, label='Cells observed before their first catalogued failure (colour and marker by forest type)'), Line2D([], [], color='#888888', marker='o', ms=4, mfc='white', mec='#888888', mew=0.8, lw=1.0, label='Never-catalogued forest')], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.09), fontsize=6.6, handlelength=2)
save(fig, 'Fig3_prefailure_stature')

# ---------------- Fig. 4: ledger ----------------
fig, ax = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True, gridspec_kw={'hspace': 0.09, 'height_ratios': [1.3, 1]})
bot = np.zeros(len(years))
for g in GS:
    v = T['loss_' + g].values; ax[0].bar(years, -v, bottom=-bot, color=COL[g], width=0.78, edgecolor='white', linewidth=0.5, label=LAB[g] + ' (loss)', zorder=3); bot += v
ax[0].errorbar(years, -T.loss, yerr=[T.loss_hi - T.loss, T.loss - T.loss_lo], fmt='none', ecolor=INK, elinewidth=0.6, capsize=1.6, zorder=4)
ax[0].bar(years, T.gain, color=COL['gain'], width=0.78, edgecolor='white', linewidth=0.5, hatch='////', label='Re-accumulation on earlier scars', zorder=3)
ax[0].errorbar(years, T.gain, yerr=[T.gain - T.gain_lo, T.gain_hi - T.gain], fmt='none', ecolor=INK, elinewidth=0.6, capsize=1.6, zorder=4)
big = T.nir_disturbance.fillna(0).values >= 0.05; ax[0].plot(years[big], -T.nir_disturbance.values[big], ls='none', marker='D', ms=4.6, mfc=COL['official'], mec='white', mew=0.6, label='Disturbance loss in the national inventory (years $\\geq$ 0.05 Mt)', zorder=5)
ax[0].axhline(0, color=INK, lw=0.6, zorder=2)
for y in (2007, 2009, 2016, 2024): ax[0].text(y, -T.loss_hi[y - 2004] - 0.22, f'{T.loss[y - 2004]:.2f}', ha='center', va='top', fontsize=6.8)
ax[0].text(2009.55, -2.535, '2.54', fontsize=6.5, va='center', color=INK)
ax[0].set_ylabel('Living-biomass stock change (Mt CO$_2$ yr$^{-1}$)'); ax[0].set_ylim(-10.6, 1.2); hgrid(ax[0]); ax[0].legend(loc='lower right', ncol=1, handlelength=1.6); panel(ax[0], 'a', x=-0.085, y=0.99)
bt = np.load(D + '/bootstrap.npz'); cum = np.cumsum(T.gain.values - T.loss.values); cb_ = np.percentile(np.cumsum(bt['gain'] - bt['loss'], 1), [2.5, 97.5], 0)
ax[1].fill_between(years, cb_[0], cb_[1], color='#D9D9D9', lw=0, zorder=2, label='95% interval'); ax[1].plot(years, cum, '-', color=INK, lw=1.3, marker='o', ms=2.6, label='Net change, central estimate', zorder=4)
ax[1].plot(years, -np.cumsum(T.loss.values), '--', color=INK, lw=1.0, label='Loss only, central estimate', zorder=3)
ax[1].plot(years[:-1], -np.cumsum(T.nir_disturbance.values[:-1]), '-', color=COL['official'], lw=1.3, marker='D', ms=2.8, mec='white', mew=0.4, label='Disturbance loss in the national inventory', zorder=4)
mn = lambda v: f'\u2212{abs(v):.1f}'
ax[1].text(2025.25, cum[-1], mn(cum[-1]), fontsize=6.8, va='center'); ax[1].text(2025.25, -T.loss.sum(), mn(T.loss.sum()), fontsize=6.8, va='center'); ax[1].text(2024.25, -np.nansum(T.nir_disturbance.values) + 0.9, mn(np.nansum(T.nir_disturbance.values)), fontsize=6.8, va='center')
ax[1].set_ylabel('Cumulative change\nsince 2004 (Mt CO$_2$)'); ax[1].set_xlabel('Year'); ax[1].set_xticks(years[::2]); ax[1].set_xlim(2003.3, 2026.3); ax[1].set_ylim(-18, 1); hgrid(ax[1]); ax[1].legend(loc='lower left', handlelength=1.8); panel(ax[1], 'b', x=-0.085, y=0.99)
save(fig, 'Fig5_ledger')

# ---------------- Fig. 5: recovery ----------------
rb = pd.read_csv(D + '/recovery_bins.csv'); sg = pd.read_csv(D + '/size_gradient.csv'); pt = pd.read_csv(D + '/persistence_2009_by_type.csv')
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), gridspec_kw={'hspace': 0.42, 'wspace': 0.32}); tt = np.linspace(0, 16, 30)
ax = axs[0, 0]
for g in GS:
    b_ = rb[(rb.tg == g) & (rb.n >= 100)]; ax.plot(b_.age, b_.x, ls='none', marker=MK[g], ms=4.5, mfc=COL[g], mec='white', mew=0.6, zorder=3)
    am = np.average(b_.age, weights=b_.n); xm = np.average(b_.x, weights=b_.n); ax.plot(tt, xm + P['k'][g] * (tt - am), '-', color=COL[g], lw=1.1, zorder=2, label=f"{LAB[g]}: {100 * P['k'][g]:.2f}% yr$^{{-1}}$")
ax.set_xlim(0, 16); ax.set_ylim(0, 1.0); ax.set_yticks(np.arange(0, 1.01, 0.2)); ax.set_xlabel('Years since failure at the image date'); ax.set_ylabel('Relative canopy height on scars\n(local reference forest = 1)'); hgrid(ax); ax.legend(loc='upper left', handlelength=1.6, fontsize=6.4, borderaxespad=0.3); panel(ax, 'a', x=-0.19)
ax = axs[0, 1]; xs = np.arange(6)
ax.errorbar(xs, 100 * sg.k, yerr=[100 * (sg.k - sg.k_lo), 100 * (sg.k_hi - sg.k)], fmt='o', color=INK, ms=4.5, mec='white', mew=0.6, capsize=2.5, elinewidth=0.8, zorder=3); ax.axhline(0, color=INK, lw=0.5)
ax.set_xticks(xs); ax.set_xticklabels(['<1', '1-3', '3-10', '10-30', '30-100', '$\\geq$100'], fontsize=6.8); ax.set_xlabel('Area of the landslide polygon (ha)'); ax.set_ylabel('Recovery rate\n(% of reference height yr$^{-1}$)'); ax.set_ylim(-1.8, 2.2); hgrid(ax); panel(ax, 'b', x=-0.19)
ax = axs[1, 0]; yy = np.arange(2017, 2025)
for g in GS:
    r = pt[pt.tg == g].iloc[0]; v = np.array([100 * r[str(y)] for y in yy]); ax.plot(yy, v, '-', color=COL[g], lw=1.2, marker=MK[g], ms=4.5, mec='white', mew=0.6, zorder=3); ax.text(2024.2, v[-1], f'{v[-1]:.0f}', fontsize=6.6, va='center')
ax.set_xlim(2016.6, 2024.9); ax.set_ylim(0, 80); ax.set_xlabel('Year of the official landslide map'); ax.set_ylabel('Scar cells of 2009 still mapped\nas landslide (%)'); hgrid(ax)
ax.legend(handles=[Line2D([], [], color=COL[g], marker=MK[g], ms=4.5, mec='white', mew=0.6, lw=1.2, label=LAB[g]) for g in GS], loc='center', bbox_to_anchor=(0.5, 0.235), ncol=2, fontsize=6.3, handlelength=1.6, columnspacing=0.8); panel(ax, 'c', x=-0.19)
ax = axs[1, 1]; w = 0.36
ax.bar(xs - w / 2, 100 * sg.bare2017, w, color='white', edgecolor=INK, linewidth=0.7, label='Map of 2017', zorder=3); ax.bar(xs + w / 2, 100 * sg.bare2024, w, color=INK, edgecolor=INK, linewidth=0.7, label='Map of 2024', zorder=3)
ax.set_xticks(xs); ax.set_xticklabels(['<1', '1-3', '3-10', '10-30', '30-100', '$\\geq$100'], fontsize=6.8); ax.set_xlabel('Area of the landslide polygon (ha)'); ax.set_ylabel('Scar cells still mapped\nas landslide (%)'); ax.set_ylim(0, 60); hgrid(ax); ax.legend(loc='upper left', ncol=2, handlelength=1.2); panel(ax, 'd', x=-0.19)
save(fig, 'Fig6_recovery')

# ---------------- Fig. 6: structural variants ----------------
V = S['variants']; order = ['central', 'unadjusted stand stock (tree cover 2000 >= 30%)', 'reference stature for failures up to 2010', 'stratum-mean stock for all cells', 'quadratic biomass-height scaling',
                           'recovery rate: within failure-year regression', 'recovery rate: within image-year and size regression', 'recovery rate: level at 8-12 yr minus baseline', 'scar-size-specific recovery', 'recovery parameters including cells with repeated failures',
                           'parameters from cells with landslide fraction >= 0.5', 'parameters from cells with landslide fraction >= 0.95', 'parameters from scar interiors (3x3 mean fraction >= 0.8)']
names = ['Central estimate', 'Stand stock without stature adjustment', 'Reference stature for failures to 2010', 'Stratum-mean stock for all cells', 'Quadratic biomass-height scaling', 'Recovery: within failure year', 'Recovery: within image year and size class', 'Recovery: level at 8-12 yr minus baseline',
         'Recovery: by scar-size class', 'Recovery: cells with repeated failures included', 'Parameters: landslide fraction $\\geq$ 0.5', 'Parameters: landslide fraction $\\geq$ 0.95', 'Parameters: scar interiors only']
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.9), sharey=True, gridspec_kw={'wspace': 0.1, 'left': 0.33}); yp = np.arange(len(order))[::-1]
for ax, key, lab, ci in [(axs[0], 'loss_total', 'Loss 2004-2025\n(Mt CO$_2$)', S['loss_total_ci']), (axs[1], 'gain_total', 'Re-accumulation 2004-2025\n(Mt CO$_2$)', S['gain_total_ci']), (axs[2], 'offset_pct', 'Re-accumulation / loss\n(%)', S['offset_pct_ci'])]:
    ax.axvspan(ci[0], ci[1], color='#E6E6E6', lw=0, zorder=1); ax.axvline(V['central'][key], color=INK, lw=0.6, zorder=2)
    v = [V[o][key] for o in order]; ax.plot(v, yp, ls='none', marker='o', ms=4.5, mfc=INK, mec='white', mew=0.6, zorder=3); ax.plot(v[0], yp[0], ls='none', marker='o', ms=6, mfc=INK, mec='white', mew=0.6, zorder=4)
    span = max(v) * 1.22
    for vv, y in zip(v, yp): ax.text(vv + 0.035 * span, y, f'{vv:.1f}', ha='left', va='center', fontsize=6.2, color=MUTE, zorder=6, bbox=dict(fc='white', ec='none', pad=0.4, alpha=0.85))
    ax.set_xlabel(lab); ax.grid(True, axis='x', color=GRID, lw=0.5); ax.set_axisbelow(True); ax.set_ylim(-0.7, len(order) - 0.3); ax.set_xlim(0, span)
axs[0].set_yticks(yp); axs[0].set_yticklabels(names, fontsize=6.8)
save(fig, 'Fig8_variants')

# ---------------- supplementary figures ----------------
# S1: image years and reference height
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={'wspace': 0.3}); ys = sorted(float(k) for k in Q['chm_year_share'])
axs[0].bar(ys, [100 * Q['chm_year_share'][str(y)] for y in ys], color='#636363', width=0.75, edgecolor='white', linewidth=0.5); axs[0].set_xlabel('Image year of the canopy height map'); axs[0].set_ylabel('Share of landslide cells (%)'); hgrid(axs[0]); panel(axs[0], 'a', x=-0.17)
ry = Q['ref_height_by_image_year']; axs[1].plot([float(k) for k in ry], list(ry.values()), '-', color=INK, lw=1.1, marker='o', ms=4, mec='white', mew=0.6); axs[1].set_ylim(0, 19); axs[1].set_xlabel('Image year of the canopy height map'); axs[1].set_ylabel('Mean canopy height of\nnever-catalogued forest (m)'); hgrid(axs[1]); panel(axs[1], 'b', x=-0.2)
save(fig, 'FigS1_image_years')
# S2: GEDI
gt = pd.read_csv(D + '/gedi_agbd_by_height_slope_BL.csv'); gsc = pd.read_csv(D + '/gedi_scar_agbd_by_age_slope.csv')
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={'wspace': 0.28}); ramp = ['#C7C7C7', '#A0A0A0', '#787878', '#4F4F4F', '#1A1A1A']; sl = ['<10°', '10-20°', '20-30°', '30-40°', '$\\geq$40°']
for s_ in range(5):
    d = gt[(gt.sb == s_) & (gt.n >= 100)]; axs[0].plot(d.h, d.agbd, '-', color=ramp[s_], lw=1.2, marker='o', ms=3.6, mec='white', mew=0.5, label=sl[s_])
axs[0].set_xlim(0, 30.5); axs[0].set_ylim(0, 310); axs[0].set_xlabel('CHMv2 mean canopy height at the footprint (m)'); axs[0].set_ylabel('GEDI L4A biomass density (Mg ha$^{-1}$)'); hgrid(axs[0]); axs[0].legend(title='Terrain slope', loc='lower right', ncol=2, fontsize=6.2, title_fontsize=6.4, handlelength=1.4); panel(axs[0], 'a', x=-0.17)
r3 = ['#BDBDBD', '#737373', '#252525']; s3 = ['<20°', '20-30°', '$\\geq$30°']; ab = ['(0, 2]', '(2, 5]', '(5, 10]', '(10, 13]', '(13, 16]', '(16, 22]']; xa = np.arange(6)
for s_ in range(3):
    d = gsc[gsc.s3 == s_].set_index('ab').reindex(ab); axs[1].plot(xa, d.agbd, '-', color=r3[s_], lw=1.2, marker='o', ms=3.8, mec='white', mew=0.5, label=s3[s_])
axs[1].set_xticks(xa); axs[1].set_xticklabels(['0-2', '2-5', '5-10', '10-13', '13-16', '16-22'], fontsize=6.8); axs[1].set_ylim(0, 310); axs[1].set_xlabel('Years since failure at the GEDI shot'); axs[1].set_ylabel('GEDI L4A biomass density on scars (Mg ha$^{-1}$)'); hgrid(axs[1]); axs[1].legend(title='Terrain slope', loc='lower right', fontsize=6.2, title_fontsize=6.4, handlelength=1.4); panel(axs[1], 'b', x=-0.17)
save(fig, 'FigS3_gedi')
# S3: persistence matrix
per = pd.read_csv(D + '/persistence_by_failure_year.csv'); M = 100 * per[[str(y) for y in range(2017, 2025)]].values
fig, ax = plt.subplots(figsize=(5.0, 3.6)); im = ax.imshow(M, cmap='Greys', vmin=0, vmax=80, aspect='auto')
ax.set_xticks(range(8)); ax.set_xticklabels(range(2017, 2025)); ax.set_yticks(range(len(per))); ax.set_yticklabels([f"{int(y)}  (n = {int(n):,})".replace(",", "\u2009") for y, n in zip(per['last'], per.n)], fontsize=6.5)
for i in range(M.shape[0]):
    for j in range(M.shape[1]): ax.text(j, i, f'{M[i, j]:.0f}', ha='center', va='center', fontsize=6.2, color='white' if M[i, j] > 45 else INK)
ax.set_xlabel('Year of the official landslide map'); ax.set_ylabel('Year of failure (cells with one catalogued failure)'); cb = fig.colorbar(im, ax=ax, pad=0.02); cb.set_label('Cells still mapped as landslide (%)'); cb.outline.set_linewidth(0.5)
for sp in ax.spines.values(): sp.set_visible(False)
save(fig, 'FigS2_persistence_matrix')
# S4: recovery rate specifications
ks = P['k_specs']; sp_ = ['elev', 'elev+size', 'elev+image year', 'elev+image year+size', 'elev+failure year', 'elev+failure year+size', 'elev+block+size']; sn = ['Elevation band', '+ size class (central)', '+ image year', '+ image year + size class', '+ failure year', '+ failure year + size class', '+ 10 km block + size class']
fig, ax = plt.subplots(figsize=(5.6, 3.0)); yp = np.arange(len(sp_))[::-1]
for g, dy in zip(GS[:3], (0.18, 0, -0.18)): ax.plot([100 * ks[s_][g] for s_ in sp_], yp + dy, ls='none', marker=MK[g], ms=4.6, mfc=COL[g], mec='white', mew=0.6, label=LAB[g])
ax.set_yticks(yp); ax.set_yticklabels(sn, fontsize=6.8); ax.set_xlim(0, 1.6); ax.set_xlabel('Gain of relative canopy height (% of reference height yr$^{-1}$)'); ax.grid(True, axis='x', color=GRID, lw=0.5); ax.set_axisbelow(True); ax.legend(loc='center right', fontsize=6.4)
save(fig, 'FigS4_rate_specifications')
print('figures written:', sorted(os.listdir(FIG)))
