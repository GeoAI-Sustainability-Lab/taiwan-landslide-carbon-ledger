# -*- coding: utf-8 -*-
"""s14: independent check of pre-failure stand condition with the 3rd National Forest Inventory photo points (1990-1994).
Each photo point (500 m staggered grid, aerial-photo interpretation) carries a canopy-density class (DENSITY 0-4) and a
growing-stock class (VOLUME 0-5: 0 non-forest, 1 0-49, 2 50-99, 3 100-199, 4 200-399, 5 >=400 m3/ha).
Points inside later catalogue scars are compared with points in never-catalogued forest of the same 10-km block and the same
2000 tree-cover class, separately for failure cohorts. Output: v2/nfi3_check.csv, v2/nfi3_check.json."""
import json, numpy as np, pandas as pd, rasterio, cv2
from rasterio.windows import Window
from config import SRC as DS, CUR, WORK
D = WORK + '/v2'
X0, Y1, RES, W, H = 151050.0, 2798910.0, 30.0, 6680, 12558
rng = np.random.default_rng(20261118)

def aligned(path, band):
    with rasterio.open(path) as ds:
        c0 = int(round((X0 - ds.bounds.left) / RES)); r0 = int(round((ds.bounds.top - Y1) / RES))
        a = ds.read(band, window=Window(c0, r0, W, H), boundless=True, fill_value=np.nan).astype(np.float32)
    a[~np.isfinite(a)] = np.nan
    return a

P = pd.read_csv(DS + '/Forestry_Reference/Taiwan_NFI3_SamplePoints_FANC_MapServer_20261111/compiled/nfi3_photo_points_layer1_core.csv')
P['yr'] = P.DATE // 10000 + 1911
P = P[(P.yr >= 1990) & (P.yr <= 1994)].copy()
P['col'] = np.floor((P.geom_x - X0) / RES).astype(int); P['row'] = np.floor((Y1 - P.geom_y) / RES).astype(int)
P = P[(P.col >= 0) & (P.col < W) & (P.row >= 0) & (P.row < H)].copy()
P['cell'] = P.row * W + P.col
VMID = {0: 0.0, 1: 25.0, 2: 75.0, 3: 150.0, 4: 300.0, 5: 500.0}
P['vol'] = P.VOLUME.map(VMID); P = P[P.vol.notna()]
P['dense'] = (P.DENSITY >= 3).astype(float); P['nonforest'] = (P.VOLUME == 0).astype(float)

TC = aligned(DS + '/Earth_Observation/Lab_platform_assets/hansen_gfc_2025_30m.tif', 1)
T = np.load(CUR + '/Forestry_Reference/Taiwan_forest_type_map_30m_TWD97/forest_type_30m.npz')['type']
P['tc2000'] = TC.ravel()[P.cell.values]; P['ftype'] = T.ravel()[P.cell.values]
P['tcb'] = np.digitize(P.tc2000.values, [25, 50, 75, 90])

C = pd.read_parquet(D + '/cells.parquet', columns=['cell', 'is_ref', 'first', 'nyears', 'a', 'a3', 'tg'])
L = C[C.is_ref == 0]
cat = np.zeros(H * W, np.uint8); cat[L.cell.values] = 1
cat = cv2.dilate(cat.reshape(H, W), np.ones((5, 5), np.uint8)).ravel()      # 2-cell buffer around any catalogued cell
forest = np.isin(T.ravel(), [1, 2, 3, 4, 5, 6])

S = P.merge(L[['cell', 'first', 'nyears', 'a', 'a3', 'tg']], on='cell', how='inner')
S = S[(S.a3 >= 0.5)].copy()
R = P[(cat[P.cell.values] == 0) & forest[P.cell.values]].copy()
for X in (S, R): X['blk'] = (X.row // 333) * 100 + (X.col // 333)
ref = R.groupby(['blk', 'tcb']).agg(rv=('vol', 'mean'), rd=('dense', 'mean'), rn=('nonforest', 'mean'), nref=('vol', 'size')).reset_index()
S = S.merge(ref[ref.nref >= 10], on=['blk', 'tcb'], how='inner')

def cohort(y):
    return '2004–2008' if y <= 2008 else ('2009' if y == 2009 else ('2010–2015' if y <= 2015 else '2016–2025'))
S['coh'] = S['first'].astype(int).map(cohort)

def stats(d):
    o = {'n_points': int(len(d)), 'n_blocks': int(d.blk.nunique()),
         'vol_scar': float(d.vol.mean()), 'vol_ref': float(d.rv.mean()), 'rel_vol': float(d.vol.mean() / d.rv.mean()),
         'dense_scar': float(d.dense.mean()), 'dense_ref': float(d.rd.mean()), 'nonforest_scar': float(d.nonforest.mean()), 'nonforest_ref': float(d.rn.mean())}
    blks = d.blk.unique(); g = {b: x for b, x in d.groupby('blk')}; bs = []
    for _ in range(500):
        pick = rng.choice(blks, len(blks)); x = pd.concat([g[b] for b in pick])
        bs.append(x.vol.mean() / x.rv.mean())
    o['rel_vol_lo'], o['rel_vol_hi'] = [float(v) for v in np.percentile(bs, [2.5, 97.5])]
    return o

rows = []
for sub, lab in ((S, 'all tree cover'), (S[S.tcb >= 3], 'tree cover 2000 >=75%')):
    for c in ['2004–2008', '2009', '2010–2015', '2016–2025']:
        d = sub[sub.coh == c]
        if len(d) >= 30: rows.append(dict(subset=lab, cohort=c, **stats(d)))
    d = sub[sub['first'] >= 2011]
    rows.append(dict(subset=lab, cohort='2011–2025', **stats(d)))
out = pd.DataFrame(rows); out.to_csv(D + '/nfi3_check.csv', index=False)

# rho used by the ledger for the same cohorts (cell-level assigned pre-failure relative biomass, tree-cover >=75%)
Pm = json.load(open(D + '/parameters.json'))
G = {'BL': 0, 'CF': 1, 'MX': 2, 'BB': 3}
def rho_assigned(d):
    v = [Pm['rho'][g][int(t)] for g, t in zip(d.tg, d.tcb) if isinstance(g, str) and g in Pm['rho']]
    return float(np.mean(v)) if v else float('nan')
extra = {}
for c in ['2004–2008', '2009', '2010–2015', '2016–2025']:
    d = S[(S.coh == c) & (S.tcb >= 3)]
    extra[c] = {'rho_assigned_tc75': rho_assigned(d)}
info = {'n_points_total': int(len(P)), 'n_points_scar': int(len(S)), 'n_points_ref': int(len(R)), 'years': sorted(P.yr.unique().tolist()),
        'rho_assigned': extra}
json.dump(info, open(D + '/nfi3_check.json', 'w'), indent=1, ensure_ascii=False)
pd.set_option('display.width', 220)
print(out.round(3).to_string()); print(json.dumps(info, ensure_ascii=False))
