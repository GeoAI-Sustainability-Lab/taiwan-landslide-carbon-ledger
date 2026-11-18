# -*- coding: utf-8 -*-
"""s02: cell-level table for the national landslide carbon ledger.
One row per 30 m cell touched by any catalogued polygon (2004-2025), plus reference forest cells.
Input locations are set in config.py; the output directory is the first argument (v1/ in the working directory).
"""
import os, sys, json, glob
import numpy as np, pandas as pd, rasterio, pyproj
from rasterio.windows import Window
from config import SRC as DS, CUR, WORK
Q1 = WORK + '/q01'
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
X0, Y1, RES, W, H = 151050.0, 2798910.0, 30.0, 6680, 12558
rng = np.random.default_rng(20261111)

T = np.load(CUR + '/Forestry_Reference/Taiwan_forest_type_map_30m_TWD97/forest_type_30m.npz')['type']
V = np.load(CUR + '/Forestry_Reference/Taiwan_4th_NFI_volume_30m_TWD97/nfi4_volume_30m.npz')['volume']
DAYS = np.load(WORK + '/s01/chmv2_acqdate_30m.npz')['days']

def aligned(path, band):
    """read a band of an EPSG:3826 30 m raster onto the project grid"""
    with rasterio.open(path) as ds:
        c0 = int(round((X0 - ds.bounds.left) / RES)); r0 = int(round((ds.bounds.top - Y1) / RES))
        a = ds.read(band, window=Window(c0, r0, W, H), boundless=True, fill_value=np.nan).astype(np.float32)
    a[~np.isfinite(a)] = np.nan
    return a
TC = aligned(DS + '/Earth_Observation/Lab_platform_assets/hansen_gfc_2025_30m.tif', 1)
ELEV = aligned(DS + '/Earth_Observation/Lab_platform_assets/dem_terrain_30m.tif', 1)
SLOPE = aligned(DS + '/Earth_Observation/Lab_platform_assets/dem_terrain_30m.tif', 2)

# ---- events per (cell, year)
z = np.load(Q1 + '/q01_record_cells.npz')
ev = pd.DataFrame({'fid': z['fid'], 'cell': z['cell'].astype(np.int64), 'frac': z['frac'], 'year': z['year'].astype(int)})
rec = pd.read_csv(Q1 + '/q01_records.csv', usecols=['fid', 'catalog_year', 'event_category', 'before_date', 'after_date'])
ad = pd.to_datetime(rec.after_date, errors='coerce'); bd = pd.to_datetime(rec.before_date, errors='coerce')
dec = lambda d: d.dt.year + (d.dt.dayofyear - 1) / 365.25
mid = (dec(ad) + dec(bd)) / 2
ok = mid.notna() & (dec(ad) >= rec.catalog_year - 0.5) & (dec(bd) <= rec.catalog_year + 1.0) & ((dec(ad) - dec(bd)) < 2.5)
rec['t_event'] = np.where(ok, mid, rec.catalog_year + 0.62)     # default: mid-August of the catalogue year (typhoon season)
rec['t_after'] = np.where(ok, dec(ad), rec.catalog_year + 1.0)
rec['t_before'] = np.where(ok, dec(bd), rec.catalog_year + 0.0)
ev = ev.merge(rec[['fid', 't_event', 't_after', 't_before', 'event_category']], on='fid')
# several polygons of one year in a cell: sum of fractions capped at 1
cy = ev.groupby(['cell', 'year']).agg(frac=('frac', 'sum'), t_event=('t_event', 'mean'), t_after=('t_after', 'max'), t_before=('t_before', 'min'),
                                      cat=('event_category', 'first')).reset_index()
cy['frac'] = cy.frac.clip(upper=1.0)
cy.to_parquet(OUT + '/cell_year_events.parquet')
g = cy.groupby('cell').agg(a=('frac', 'max'), first=('year', 'min'), last=('year', 'max'), nyears=('year', 'nunique'),
                           t_first_before=('t_before', 'min'), t_last=('t_event', 'max'), t_last_after=('t_after', 'max')).reset_index()

def neigh_type(cells):
    r = cells // W; c = cells % W; out = T.ravel()[cells].copy()
    for i in np.nonzero((out == 0) | (out > 6))[0]:
        w = T[max(r[i] - 2, 0):r[i] + 3, max(c[i] - 2, 0):c[i] + 3].ravel(); w = w[(w > 0) & (w < 7)]
        if len(w): out[i] = np.bincount(w).argmax()
    return out
ls_mask = np.zeros(H * W, bool); ls_mask[g.cell.values] = True
cand = np.nonzero((T.ravel() >= 1) & (T.ravel() <= 6) & (~ls_mask))[0]
ref = pd.DataFrame({'cell': rng.choice(cand, 500000, replace=False)}); ref['is_ref'] = True
g['is_ref'] = False
A = pd.concat([g, ref], ignore_index=True)
cells = A.cell.values
A['type_cell'] = T.ravel()[cells]; A['type_nb'] = neigh_type(cells); A['vol'] = V.ravel()[cells]
A['tc2000'] = TC.ravel()[cells]; A['elev'] = ELEV.ravel()[cells]; A['slope'] = SLOPE.ravel()[cells]
A['t_chm'] = np.where(DAYS.ravel()[cells] > 0, 2000 + DAYS.ravel()[cells] / 365.25, np.nan)
r = cells // W; c = cells % W
A['x'] = X0 + (c + 0.5) * RES; A['y'] = Y1 - (r + 0.5) * RES
lon, lat = pyproj.Transformer.from_crs(3826, 4326, always_xy=True).transform(A.x.values, A.y.values); A['lon'] = lon; A['lat'] = lat

def sample(paths, src_crs, win, nodata, conv=None):
    px, py = pyproj.Transformer.from_crs(3826, src_crs, always_xy=True).transform(A.x.values, A.y.values)
    val = np.full(len(A), np.nan, np.float32)
    for p in paths:
        with rasterio.open(p) as ds:
            col, row = (~ds.transform) * (np.asarray(px), np.asarray(py))
            col = np.floor(col).astype(np.int64); row = np.floor(row).astype(np.int64)
            idx = np.nonzero((col >= win) & (col < ds.width - win) & (row >= win) & (row < ds.height - win))[0]
            if len(idx) == 0: continue
            STR = 2048
            for r0 in range(int(row[idx].min()) // STR * STR, int(row[idx].max()) + 1, STR):
                sel = idx[(row[idx] >= r0) & (row[idx] < r0 + STR)]
                if len(sel) == 0: continue
                c0 = int(col[sel].min()) - win; c1 = int(col[sel].max()) + win + 1
                a = ds.read(1, window=Window(c0, r0 - win, c1 - c0, STR + 2 * win), boundless=True, fill_value=nodata).astype(np.float32)
                a[a == nodata] = np.nan
                if conv: a = conv(a)
                rr = row[sel] - (r0 - win); cc = col[sel] - c0
                acc = np.zeros(len(sel), np.float32); cnt = np.zeros(len(sel), np.float32)
                for dr in range(-win, win + 1):
                    for dc in range(-win, win + 1):
                        v = a[rr + dr, cc + dc]; m = np.isfinite(v); acc[m] += v[m]; cnt[m] += 1
                val[sel] = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    return val
chm = sorted(glob.glob(DS + '/Terrain_Canopy/CHMv2_10m_mean*/mean/*.tif'))
print('CHM tiles', [os.path.basename(p) for p in chm])
A['chm'] = sample(chm, 3857, 1, 255)
# ALOS-2 PALSAR-2 2020 mosaic, HV: gamma0 (linear power) averaged over 3x3 (about 75 m), stored in dB
hv = sorted(glob.glob(DS + '/Earth_Observation/ALOS2_PALSAR2_Mosaic_CentralTaiwan_2018_2020_20261107/2020/*_HV_*.tif'))
lin = sample(hv, 4326, 1, 0, conv=lambda d: np.where(d > 1, (d.astype(np.float64) ** 2) * 10 ** (-8.3), np.nan).astype(np.float32))
A['hv2020_db'] = 10 * np.log10(lin)
A.to_parquet(OUT + '/cells.parquet')
# national strata of pre-event stand volume: forest type x elevation band, cells with a volume estimate
band = np.digitize(ELEV, [500, 1500, 2500])
rows = []
for t in range(1, 7):
    for b in range(4):
        m = (T == t) & (band == b) & np.isfinite(ELEV); vv = V[m]
        rows.append({'type': t, 'band': b, 'ha': float(m.sum()) * 0.09, 'ha_vol': float((vv > 0).sum()) * 0.09,
                     'vol_mean': float(vv[vv > 0].mean()) if (vv > 0).any() else np.nan, 'vol_sd': float(vv[vv > 0].std()) if (vv > 0).any() else np.nan})
pd.DataFrame(rows).to_csv(OUT + '/nfi4_strata_type_elev.csv', index=False)
print(len(g), 'landslide cells;', int(A.chm.notna().sum()), 'with CHM;', int(A.t_chm.notna().sum()), 'with CHM date;', int(A.hv2020_db.notna().sum()), 'with HV;',
      'tc2000 valid', int(A.tc2000.notna().sum()))
print(pd.DataFrame(rows).round(1).to_string())
