# -*- coding: utf-8 -*-
"""s22: case-site figures (Fig. 4 Heping, Hualien; Fig. 7 Laonong River, Kaohsiung) and their statistics (v2/case_stats.json).
All panels are computed from data: Landsat C2 L2 / Sentinel-2 L2A chips (Planetary Computer), CHMv2 10 m mean height,
the event catalogue (cell fractions), the official annual landslide maps (cell fractions), and the cell table of the ledger."""
import json, glob, sys, numpy as np, pandas as pd, pyproj, shapefile
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cases')); from common import chm10
from config import SRC
INK = '#222222'; MUTE = '#666666'; GRID = '#E3E3E3'; REF = '#888888'; LARGE = '#F0E442'; SMALL = '#882255'
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.6, 'axes.edgecolor': INK, 'axes.labelcolor': INK, 'text.color': INK,
                     'xtick.color': INK, 'ytick.color': INK, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'legend.fontsize': 6.8,
                     'axes.grid': False, 'savefig.dpi': 300})
W = 6680; X0, Y1 = 151050.0, 2798910.0; FIG = 'figures'; TOWGS = pyproj.Transformer.from_crs(3826, 4326, always_xy=True)
RATIO = {1: 1.9622, 2: 1.0999, 4: 1.5310, 5: 2.5256, 6: 2.1150}
def nfmt(n): return f'{int(n):,}'.replace(',', '\u2009')
def fpanel(fig, xin, yin, s, FW, FH): fig.text(xin / FW, yin / FH, s, fontsize=10, fontweight='bold', va='bottom', ha='left')
def panel(ax, s, x=-0.02, y=1.02): ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
def hgrid(ax, axis='y'): ax.grid(True, axis=axis, color=GRID, lw=0.5, ls='-'); ax.set_axisbelow(True)
def rgb(b, hi=0.14, lo=0.0, gam=0.85):
    x = np.stack([b[2], b[1], b[0]], -1); x = np.nan_to_num(x, nan=0.0); return np.clip((x - lo) / (hi - lo), 0, 1) ** gam
def up(m, f): return np.kron(m.astype(float), np.ones((f, f)))
def halo(c, w=2.2): return [pe.Stroke(linewidth=w, foreground=c), pe.Normal()]
def outline(ax, m, f, col, hc, lw=0.8):
    cs = ax.contour(up(m, f), [0.5], colors=[col], linewidths=lw, origin='upper', extent=(-0.5, m.shape[1] * f - 0.5, m.shape[0] * f - 0.5, -0.5))
    for c in cs.collections if hasattr(cs, 'collections') else [cs]: c.set_path_effects(halo(hc, lw + 1.3))
def graticule(ax, R, f, step, left=True, bottom=True, xevery=1):
    r0, r1, c0, c1 = R; h, w = (r1 - r0) * f, (c1 - c0) * f; res = 30 / f
    xs = X0 + c0 * 30 + (np.arange(w) + 0.5) * res; ys = Y1 - r0 * 30 - (np.arange(h) + 0.5) * res
    XX, YY = np.meshgrid(xs[::4], ys[::4]); lo, la = TOWGS.transform(XX, YY)
    lv = np.arange(np.ceil(lo.min() / step) * step, lo.max(), step); av = np.arange(np.ceil(la.min() / step) * step, la.max(), step)
    ext = (-0.5, w - 0.5, h - 0.5, -0.5)
    ax.contour(lo, lv, colors='white', linewidths=0.35, alpha=0.75, extent=(ext[0], ext[1], ext[2], ext[3]), origin='upper')
    ax.contour(la, av, colors='white', linewidths=0.35, alpha=0.75, extent=(ext[0], ext[1], ext[2], ext[3]), origin='upper')
    bx = np.interp(lv, lo[-1], np.arange(lo.shape[1]) * 4); ly = np.interp(av, la[::-1, 0], (np.arange(la.shape[0]) * 4)[::-1])
    ax.set_xticks(bx); ax.set_yticks(ly)
    ax.set_xticklabels([f'{v:.2f}°E' if i % xevery == 0 else '' for i, v in enumerate(lv)] if bottom else []); ax.set_yticklabels([f'{v:.2f}°N' for v in av] if left else [], rotation=90, va='center')
    ax.tick_params(length=2, labelsize=5.8, pad=1.5); ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    for s in ax.spines.values(): s.set_visible(True)
def scalebar(ax, f, km=1.0, x=0.05, y=0.06):
    px = km * 1000 / (30 / f); x0 = ax.get_xlim()[0] + x * (ax.get_xlim()[1] - ax.get_xlim()[0]); y0 = ax.get_ylim()[0] - y * (ax.get_ylim()[0] - ax.get_ylim()[1])
    ax.plot([x0, x0 + px], [y0, y0], color='white', lw=2.2, solid_capstyle='butt', path_effects=halo(INK, 3.4)); ax.text(x0 + px / 2, y0 - 0.035 * (ax.get_ylim()[0] - ax.get_ylim()[1]), f'{km:g} km', color='white', ha='center', va='bottom', fontsize=6.2, path_effects=halo(INK, 1.6))
def datelab(ax, s): ax.text(0.03, 0.97, s, transform=ax.transAxes, color='white', fontsize=6.4, va='top', path_effects=halo(INK, 1.8))
sf = shapefile.Reader(SRC + '/Reference_Boundaries/NaturalEarth_admin_10m/ne_10m_admin_0_countries.shp', encoding='latin1')
fld = [f[0] for f in sf.fields[1:]]; ia = fld.index('ADM0_A3'); TW = [s for s, r in zip(sf.shapes(), sf.records()) if r[ia] == 'TWN'][0]
def locator(ax, sites):
    pts = np.array(TW.points); parts = list(TW.parts) + [len(pts)]
    for k in range(len(parts) - 1):
        p = pts[parts[k]:parts[k + 1]]
        if p[:, 0].min() > 119.9 and p[:, 1].min() > 21.8: ax.plot(p[:, 0], p[:, 1], color=INK, lw=0.5)
    ax.set_xlim(119.95, 122.05); ax.set_ylim(21.85, 25.35); ax.set_aspect(1 / np.cos(np.radians(23.6)))
    ax.set_xticks([120, 121, 122]); ax.set_yticks([22, 23, 24, 25]); ax.set_xticklabels(['120°E', '121°E', '122°E']); ax.set_yticklabels(['22°N', '23°N', '24°N', '25°N'])
    ax.grid(True, color='#BBBBBB', lw=0.4); ax.tick_params(labelsize=5.8, length=2)
    for s in ax.spines.values(): s.set_visible(True)
    for lab, lo, la, hl in sites:
        ax.plot(lo, la, 's', ms=5.5 if hl else 4, mfc=INK if hl else 'white', mec=INK, mew=0.8)
        ax.text(lo + 0.08, la, lab, fontsize=6.3, va='center', fontweight='bold' if hl else 'normal')

C = pd.read_parquet('v2/cells.parquet'); C['row'] = C.cell // W; C['col'] = C.cell % W
E = pd.read_parquet('v1/cell_year_events.parquet'); E['row'] = E.cell // W; E['col'] = E.cell % W
OUT = {}
# =========================================================== site A: Laonong River, Taoyuan District, Kaohsiung (Morakot)
RA = (7500, 7800, 2300, 2700); r0, r1, c0, c1 = RA; shp = (r1 - r0, c1 - c0)
c = C[(C.row >= r0) & (C.row < r1) & (C.col >= c0) & (C.col < c1)]
s9 = c[(c.is_ref == 0) & (c['first'] == 2009) & (c.a >= 0.8)]
G = {'large': s9[s9.poly_ha >= 100], 'small': s9[(s9.poly_ha < 3) & (s9.nyears == 1) & (s9.tc2000 >= 75)], 'ref': c[c.is_ref == 1]}
M = {}
for k in ('large', 'small'):
    m = np.zeros(shp, bool); m[G[k].row - r0, G[k].col - c0] = True; M[k] = m
chA = np.load('cases/A_chm10.npy')
ts = pd.concat([pd.read_csv(f) for f in sorted(glob.glob('cases/A_ndvi_20*_20*.csv'))]); ts = ts[ts.ref > 0.5].copy()
for k in ('large', 'small'): ts[k + '_r'] = ts[k] / ts.ref
ts['t'] = pd.to_datetime(ts.date); ts['yr'] = ts.t.dt.year; ts['dec'] = ts.yr + (ts.t.dt.dayofyear - 0.5) / 365.25
ann = ts.groupby('yr')[['large_r', 'small_r']].median(); ts.to_csv('cases/A_ndvi_all_scenes.csv', index=False)
lon0, lat0 = TOWGS.transform(X0 + (c0 + c1) / 2 * 30, Y1 - (r0 + r1) / 2 * 30)
A = {'lon': lon0, 'lat': lat0, 'n_scenes': int(len(ts)), 'n_years': int(ts.yr.nunique())}
for k in ('large', 'small', 'ref'):
    v = G[k]; A[k] = {'n': int(len(v)), 'ha': float(v.a.sum() * 0.09), 'chm_mean': float(v.chm.mean()), 'chm_median': float(v.chm.median()), 'lt2': float((v.chm < 2).mean()), 't_chm_median': float(v.t_chm.median()),
                      'map': {str(y): float((v[f's{y}'] >= 0.5).mean()) for y in range(2017, 2025)}, 'poly_min': float(v.poly_ha.min()) if k != 'ref' else None, 'poly_max': float(v.poly_ha.max()) if k != 'ref' else None,
                      'elev_median': float(v.elev.median())}
for k in ('large', 'small'):
    A[k]['ndvi_r'] = {str(y): float(ann.loc[y, k + '_r']) for y in ann.index}
v = G['small'][np.floor(G['small'].t_chm).isin([2014, 2015])]; A['small_1415'] = {'n': int(len(v)), 'chm_mean': float(v.chm.mean()), 'lt2': float((v.chm < 2).mean())}
A['big_poly_first_years'] = sorted(set(int(y) for y in E[E.cell.isin(G['large'].cell)].year))
A['large_share_catalogued_before_2009'] = float(E[E.cell.isin(G['large'].cell) & (E.year < 2009)].cell.nunique() / len(G['large']))
A['large_recat_2016_2025'] = float(E[E.cell.isin(G['large'].cell) & (E.year >= 2016)].cell.nunique() / len(G['large']))
A['large_recat_2024_2025'] = float(E[E.cell.isin(G['large'].cell) & (E.year >= 2024)].cell.nunique() / len(G['large']))
A['small_recat_2016_2025'] = float(E[E.cell.isin(G['small'].cell) & (E.year >= 2016)].cell.nunique() / len(G['small']))
OUT['A'] = A
L5 = {k: np.load(f'cases/A_{k}.npz')['b'] for k in ['LT05_L2SP_117044_20080824_02_T1', 'LT05_L2SP_117044_20090912_02_T1']}
S2 = {k: np.load(glob.glob(f'cases/A_{k}*.npz')[0])['b'] for k in ['S2B_MSIL2A_20171117', 'S2B_MSIL2A_20240105']}
imgs = [(rgb(L5['LT05_L2SP_117044_20080824_02_T1']), 1, '24 Aug 2008, Landsat 5'), (rgb(L5['LT05_L2SP_117044_20090912_02_T1']), 1, '12 Sep 2009, Landsat 5'),
        (rgb(S2['S2B_MSIL2A_20171117']), 3, '17 Nov 2017, Sentinel-2'), (rgb(S2['S2B_MSIL2A_20240105']), 3, '5 Jan 2024, Sentinel-2')]

FW, FH = 7.2, 3.86; Lm, Rm, GAP = 0.70, 0.06, 0.40; WCOL = (FW - Lm - Rm - 3 * GAP) / 4
def cx(j, span=1): return Lm + j * (WCOL + GAP), span * WCOL + (span - 1) * GAP          # left edge and width (inches) of column j
fig = plt.figure(figsize=(FW, FH))
IH = WCOL * shp[0] / shp[1]; Y1T = FH - 0.24 - IH; CH = 1.38; Y2 = 0.80
LX = [cx(0)[0] - 0.62] + [cx(j)[0] - 0.02 for j in (1, 2, 3)]               # one letter position per column
for j, (im, f, lab) in enumerate(imgs):
    x0, w = cx(j); ax = fig.add_axes([x0 / FW, Y1T / FH, w / FW, IH / FH]); ax.imshow(im, interpolation='nearest', aspect='auto')
    outline(ax, M['large'], f, LARGE, INK, 0.7); outline(ax, M['small'], f, SMALL, 'white', 0.6)
    graticule(ax, RA, f, 0.04, left=(j == 0)); datelab(ax, lab); fpanel(fig, LX[j], Y1T + IH + 0.03, 'abcd'[j], FW, FH)
    if j == 0: scalebar(ax, f)
# e: relative NDVI (columns 1-2)
x0, w = cx(0, 2); ax = fig.add_axes([x0 / FW, Y2 / FH, w / FW, CH / FH]); ax.axhline(1, color=REF, lw=1.0, zorder=2)
ax.axvline(2009.6, color=MUTE, lw=0.6, ls='--', zorder=1); ax.text(2009.35, 0.05, 'Typhoon Morakot', fontsize=6, color=MUTE, va='bottom', ha='right')
for k, col, hc in (('large', LARGE, INK), ('small', SMALL, None)):
    ax.plot(ts.dec, ts[k + '_r'], 'o', ms=1.8, color=col, mec=INK if k == 'large' else col, mew=0.25, alpha=0.75, zorder=3)
    ln, = ax.plot(ann.index + 0.5, ann[k + '_r'], '-', color=col, lw=1.5, zorder=4)
    if hc: ln.set_path_effects(halo(hc, 2.4))
ax.set_xlim(2005, 2026); ax.set_ylim(0, 1.1); ax.set_xticks(range(2005, 2026, 5)); ax.set_ylabel('NDVI ratio to\nnever-catalogued forest'); ax.set_xlabel('Year'); hgrid(ax)
ax.text(2025.8, 0.04, 'points: single scenes · lines: annual medians', fontsize=5.6, color=MUTE, ha='right'); fpanel(fig, LX[0], Y2 + CH + 0.04, 'e', FW, FH)
# f: canopy height by group (column 3)
x0, w = cx(2); ax = fig.add_axes([x0 / FW, Y2 / FH, w / FW, CH / FH]); data = [G[k].chm.dropna().values for k in ('large', 'small', 'ref')]
bp = ax.boxplot(data, positions=[0, 1, 2], widths=0.55, whis=(5, 95), patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.0), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
for bx_, col in zip(bp['boxes'], (LARGE, SMALL, REF)): bx_.set_facecolor(col); bx_.set_edgecolor(INK)
for i, k in enumerate(('large', 'small', 'ref')): ax.text(i, 33.0, f"{A[k]['chm_mean']:.1f}", ha='center', fontsize=6.2)
ax.text(1, 36.2, 'mean (m)', fontsize=5.8, color=MUTE, ha='center')
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['≥100 ha', '<3 ha', 'Forest'], fontsize=6.2); ax.set_ylim(0, 38); ax.set_yticks(range(0, 31, 10)); ax.set_ylabel('Canopy height (m)'); hgrid(ax)
ax.text(0.5, -0.2, 'CHMv2, images 2014–2017', transform=ax.transAxes, ha='center', va='top', fontsize=6.0, color=MUTE); fpanel(fig, LX[2], Y2 + CH + 0.04, 'f', FW, FH)
# g: official annual maps (column 4)
x0, w = cx(3); ax = fig.add_axes([x0 / FW, Y2 / FH, w / FW, CH / FH]); yrs = np.arange(2017, 2025)
for k, col, mk in (('large', LARGE, 'o'), ('small', SMALL, 's'), ('ref', REF, '^')):
    y = [100 * A[k]['map'][str(v)] for v in yrs]
    ln, = ax.plot(yrs, y, '-', color=col, lw=1.4, marker=mk, ms=3.4, mec=INK if k == 'large' else col, mew=0.5, mfc=col if k != 'ref' else 'white')
    if k == 'large': ln.set_path_effects(halo(INK, 2.3))
ax.set_ylim(0, 60); ax.set_xticks([2017, 2020, 2024]); ax.set_ylabel('Mapped as landslide (%)'); ax.set_xlabel('Official map year'); hgrid(ax); fpanel(fig, LX[3], Y2 + CH + 0.04, 'g', FW, FH)
fig.legend(handles=[Line2D([], [], color=LARGE, lw=1.5, marker='o', ms=3.4, mec=INK, path_effects=halo(INK, 2.4), label='Failed in 2009, polygon ≥100 ha'),
                    Line2D([], [], color=SMALL, lw=1.5, marker='s', ms=3.4, path_effects=halo('white', 2.4), label='Failed only in 2009, polygon <3 ha, tree cover 2000 ≥75%'),
                    Line2D([], [], color=REF, lw=1.5, marker='^', ms=3.4, mfc='white', label='Never-catalogued forest, sample cells')],
           loc='lower center', ncol=3, fontsize=6.2, bbox_to_anchor=(0.5, 0.0), handlelength=2.2, columnspacing=1.6, frameon=False)
fig.savefig(f'{FIG}/Fig7_case_laonong.png'); fig.savefig(f'{FIG}/Fig7_case_laonong.pdf'); plt.close(fig)

# =========================================================== site B: Heping Village, Xiulin Township, Hualien (2024 Hualien earthquake)
C24 = '#332288'   # 2024 catalogue cells (used only for this group)
RB = (3650, 3900, 5380, 5640); r0, r1, c0, c1 = RB; shp = (r1 - r0, c1 - c0)
c = C[(C.row >= r0) & (C.row < r1) & (C.col >= c0) & (C.col < c1)]
e = E[(E.row >= r0) & (E.row < r1) & (E.col >= c0) & (E.col < c1)]
F24 = np.zeros(shp); g = e[e.year == 2024]; F24[g.row - r0, g.col - c0] = g.frac
s = c[(c.is_ref == 0) & (c['first'] == 2024) & (c.a >= 0.5)].copy()
Lg = pd.read_parquet('v2/cells_ledger.parquet')[['cell', 'S_pre', 'loss_t', 'tcb']]; s = s.merge(Lg, on='cell')
s['S_stand'] = [v * RATIO[int(t)] if (v > 0 and int(t) in RATIO) else np.nan for v, t in zip(s.vol, s.type_cell)]
ref = c[c.is_ref == 1]; chB = np.load('cases/B_chm10.npy')
pre = s[s.t_chm < s.t_first_before]; mS = s.S_stand.notna(); Aa = s.a * 0.09
lonB, latB = TOWGS.transform(X0 + (c0 + c1) / 2 * 30, Y1 - (r0 + r1) / 2 * 30)
B = {'lon': lonB, 'lat': latB, 'n': int(len(s)), 'ha': float(Aa.sum()), 'pre_share': float(len(pre) / len(s)), 't_chm_median': float(pre.t_chm.median()),
     't_chm_min': float(pre.t_chm.min()), 't_chm_max': float(pre.t_chm.max()), 'ref_t_chm_min': float(ref.t_chm.min()), 'ref_t_chm_max': float(ref.t_chm.max()),
     'chm_mean': float(pre.chm.mean()), 'lt2': float((pre.chm < 2).mean()), 'lt5': float((pre.chm < 5).mean()), 'ref_n': int(len(ref)), 'ref_chm': float(ref.chm.mean()), 'ref_lt2': float((ref.chm < 2).mean()),
     'slope_median': float(s.slope.median()), 'ref_slope_median': float(ref.slope.median()),
     'map2023': float((s.s2023 >= 0.5).mean()), 'tc_lt25': float((s.tc2000 < 25).mean()),
     'chm_by_tcb': {int(k): [int(v.size), float(v.mean())] for k, v in pre.groupby('tcb').chm},
     'stock_stand_kt': float((Aa[mS] * s.S_stand[mS]).sum() / 1e3), 'stock_adj_kt': float((Aa[mS] * s.S_pre[mS]).sum() / 1e3), 'n_stand': int(mS.sum()),
     'types': {k: int(v) for k, v in s.tg.value_counts().items()}}
B['stock_ratio'] = B['stock_adj_kt'] / B['stock_stand_kt']; OUT['B'] = B
pB = np.load(glob.glob('cases/B_S2B_MSIL2A_20240214*.npz')[0])['b']; qB = np.load(glob.glob('cases/B_S2C_MSIL2A_20250325*.npz')[0])['b']
M24 = F24 >= 0.5
FW, FH = 7.2, 4.95; Lm, Rm, GAP = 0.56, 0.06, 0.62; WCOL = (FW - Lm - Rm - 2 * GAP) / 3
def cx(j): return Lm + j * (WCOL + GAP), WCOL
fig = plt.figure(figsize=(FW, FH)); IH = WCOL * shp[0] / shp[1]; Y1T = FH - 0.24 - IH; CH = 1.42; Y2 = 0.74
LX = [cx(0)[0] - 0.50] + [cx(j)[0] - 0.02 for j in (1, 2)]
yl = f"{int(np.floor(B['t_chm_min']))}–{int(np.floor(B['t_chm_max']))}"
for j, (im, lab) in enumerate([(rgb(pB), '14 Feb 2024, Sentinel-2'), (rgb(qB), '25 Mar 2025, Sentinel-2'), (None, None)]):
    x0, w = cx(j); ax = fig.add_axes([x0 / FW, Y1T / FH, w / FW, IH / FH])
    if im is not None: ax.imshow(im, interpolation='nearest', aspect='auto'); datelab(ax, lab)
    else:
        cm = plt.get_cmap('Greens').copy(); cm.set_bad('#D9D9D9'); hm = ax.imshow(np.clip(chB, 0, 30), cmap=cm, vmin=0, vmax=30, interpolation='nearest', aspect='auto'); datelab(ax, f'CHMv2, images {yl}')
        cax = ax.inset_axes([0.0, -0.20, 1.0, 0.045]); cb = fig.colorbar(hm, cax=cax, orientation='horizontal'); cb.set_label('Canopy height before the earthquake (m)', fontsize=6.2, labelpad=1); cb.ax.tick_params(labelsize=5.8, length=2); cb.outline.set_linewidth(0.4)
    if j != 1: outline(ax, M24, 3, C24, 'white', 0.55)
    graticule(ax, RB, 3, 0.02, left=(j == 0), xevery=2); fpanel(fig, LX[j], Y1T + IH + 0.03, 'abc'[j], FW, FH)
    if j == 0: scalebar(ax, 3)
# d: ECDF of pre-earthquake canopy height
x0, w = cx(0); ax = fig.add_axes([x0 / FW, Y2 / FH, w / FW, CH / FH])
for v, col, lab in ((pre.chm.dropna().values, C24, f'2024 catalogue cells, n = {nfmt(len(pre))}'), (ref.chm.dropna().values, REF, f'Never-catalogued forest, n = {nfmt(len(ref))}')):
    x = np.sort(v); ax.step(x, np.arange(1, x.size + 1) / x.size * 100, where='post', color=col, lw=1.4, label=lab)
    ax.axvline(v.mean(), color=col, lw=0.7, ls='--')
ax.set_xlim(0, 30); ax.set_ylim(0, 100); ax.set_xlabel('Canopy height before the earthquake (m)'); ax.set_ylabel('Cumulative share of cells (%)'); hgrid(ax); hgrid(ax, 'x')
fpanel(fig, LX[0], Y2 + CH + 0.04, 'd', FW, FH)
# e: height by tree cover 2000 class
x0, w = cx(1); ax = fig.add_axes([x0 / FW, Y2 / FH, w / FW, CH / FH]); kk = sorted(B['chm_by_tcb']); lab5 = ['0–25', '25–50', '50–75', '75–90', '>90']
ax.bar(range(len(kk)), [B['chm_by_tcb'][k][1] for k in kk], 0.62, color=C24, edgecolor='white', linewidth=0.5, zorder=3)
ax.axhline(B['ref_chm'], color=REF, lw=1.2, zorder=4); ax.text(-0.45, B['ref_chm'] + 0.4, 'Never-catalogued forest', color=MUTE, fontsize=5.9, ha='left', va='bottom')
for i, k in enumerate(kk): ax.text(i, 0.4, nfmt(B['chm_by_tcb'][k][0]), ha='center', fontsize=5.6, color='white', rotation=90, va='bottom')
ax.set_xticks(range(len(kk))); ax.set_xticklabels([lab5[k] for k in kk], fontsize=6.2); ax.set_xlabel('Tree cover in 2000 (%)'); ax.set_ylabel('Mean canopy height (m)'); ax.set_ylim(0, 17); hgrid(ax); fpanel(fig, LX[1], Y2 + CH + 0.04, 'e', FW, FH)
# f: stock removed, stand average versus pre-failure stature
x0, w = cx(2); ax = fig.add_axes([x0 / FW, Y2 / FH, w / FW, CH / FH]); vals = [B['stock_stand_kt'], B['stock_adj_kt']]
ax.bar([0], vals[:1], 0.55, color='white', edgecolor=C24, linewidth=0.8, hatch='////', zorder=3); ax.bar([1], vals[1:], 0.55, color=C24, edgecolor=C24, linewidth=0.8, zorder=3)
for i, v in enumerate(vals): ax.text(i, v + 3, f'{v:.1f}', ha='center', fontsize=6.4)
ax.set_xticks([0, 1]); ax.set_xticklabels(['Stand-average\nstock', 'Pre-failure\nstature'], fontsize=6.2); ax.set_ylabel('Stock removed (kt CO$_2$)'); ax.set_ylim(0, 175); hgrid(ax); fpanel(fig, LX[2], Y2 + CH + 0.04, 'f', FW, FH)
fig.legend(handles=[Line2D([], [], color=C24, lw=1.5, label=f'2024 catalogue cells, n = {nfmt(len(pre))}'), Line2D([], [], color=REF, lw=1.5, label=f'Never-catalogued forest, n = {nfmt(len(ref))}'),
                    Line2D([], [], color=MUTE, lw=0.8, ls='--', label='Group mean')], loc='lower center', ncol=3, fontsize=6.2, bbox_to_anchor=(0.5, 0.0), handlelength=2.2, columnspacing=1.8, frameon=False)
fig.savefig(f'{FIG}/Fig4_case_hualien.png'); fig.savefig(f'{FIG}/Fig4_case_hualien.pdf'); plt.close(fig)
json.dump(OUT, open('v2/case_stats.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(OUT['B'], ensure_ascii=False)[:1500])
