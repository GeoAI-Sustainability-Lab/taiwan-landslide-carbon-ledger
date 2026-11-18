# -*- coding: utf-8 -*-
"""s06: attach site attributes, CHMv2 height and landslide history to GEDI L4A footprints (quality-filtered by the service)."""
import os, glob, numpy as np, pandas as pd, rasterio, pyproj
from rasterio.windows import Window
from config import SRC as DS, CUR, WORK
X0, Y1, RES, W, H = 151050.0, 2798910.0, 30.0, 6680, 12558
z = np.load(DS + '/Earth_Observation/GEDI_L4A_V2_1_SlideRule_Taiwan_20261111/derived/gedi_l4a_taiwan_shots.npz')
G = pd.DataFrame({k: z[k] for k in z.files})
G['t'] = 1970 + G.time.values / 1e9 / 86400 / 365.25
x, y = pyproj.Transformer.from_crs(4326, 3826, always_xy=True).transform(G.longitude.values, G.latitude.values)
G['x'] = x; G['y'] = y; col = np.floor((x - X0) / RES).astype(np.int64); row = np.floor((Y1 - y) / RES).astype(np.int64)
ok = (col >= 0) & (col < W) & (row >= 0) & (row < H); G = G[ok].reset_index(drop=True); col = col[ok]; row = row[ok]; cell = row * W + col; G['cell'] = cell
zt = np.load(CUR + '/Forestry_Reference/Taiwan_forest_type_map_30m_TWD97/forest_type_30m.npz'); zv = np.load(CUR + '/Forestry_Reference/Taiwan_4th_NFI_volume_30m_TWD97/nfi4_volume_30m.npz')
G['type'] = zt['type'].ravel()[cell]; G['func'] = zt['func'].ravel()[cell]; G['vol'] = zv['volume'].ravel()[cell]
def aligned(path, band):
    with rasterio.open(path) as ds:
        c0 = int(round((X0 - ds.bounds.left) / RES)); r0 = int(round((ds.bounds.top - Y1) / RES))
        return ds.read(band, window=Window(c0, r0, W, H), boundless=True, fill_value=np.nan).astype(np.float32)
dem = DS + '/Earth_Observation/Lab_platform_assets/dem_terrain_30m.tif'
G['elev'] = aligned(dem, 1).ravel()[cell]; G['slope'] = aligned(dem, 2).ravel()[cell]
G['tc2000'] = aligned(DS + '/Earth_Observation/Lab_platform_assets/hansen_gfc_2025_30m.tif', 1).ravel()[cell]
LS = np.load('v2/ls_frac_30m.npy'); P = np.pad(LS, 1); S3 = sum(P[1 + dr:1 + dr + H, 1 + dc:1 + dc + W] for dr in (-1, 0, 1) for dc in (-1, 0, 1)) / 9.0
P2 = np.pad(LS > 0, 2); ANY5 = np.zeros((H, W), bool)
for dr in range(-2, 3):
    for dc in range(-2, 3): ANY5 |= P2[2 + dr:2 + dr + H, 2 + dc:2 + dc + W]
G['a'] = LS.ravel()[cell]; G['a3'] = S3.ravel()[cell].astype(np.float32); G['ls5'] = ANY5.ravel()[cell]
DAYS = np.load(WORK + '/s01/chmv2_acqdate_30m.npz')['days']
d = DAYS.ravel()[cell]; G['t_chm'] = np.where(d > 0, 2000 + d / 365.25, np.nan)
# CHMv2 10 m mean at the footprint: mean of the 3x3 block of 10 m pixels centred on the shot (about the 25 m footprint)
px, py = pyproj.Transformer.from_crs(4326, 3857, always_xy=True).transform(G.longitude.values, G.latitude.values)
val = np.full(len(G), np.nan, np.float32)
for p in sorted(glob.glob(DS + '/Terrain_Canopy/CHMv2_10m_mean*/mean/*.tif')):
    with rasterio.open(p) as ds:
        c_, r_ = (~ds.transform) * (px, py); c_ = np.floor(c_).astype(np.int64); r_ = np.floor(r_).astype(np.int64)
        idx = np.nonzero((c_ >= 1) & (c_ < ds.width - 1) & (r_ >= 1) & (r_ < ds.height - 1))[0]
        if len(idx) == 0: continue
        STR = 2048
        for r0 in range(int(r_[idx].min()) // STR * STR, int(r_[idx].max()) + 1, STR):
            sel = idx[(r_[idx] >= r0) & (r_[idx] < r0 + STR)]
            if len(sel) == 0: continue
            c0 = int(c_[sel].min()) - 1; c1 = int(c_[sel].max()) + 2
            a = ds.read(1, window=Window(c0, r0 - 1, c1 - c0, STR + 2), boundless=True, fill_value=255).astype(np.float32); a[a == 255] = np.nan
            rr = r_[sel] - (r0 - 1); cc = c_[sel] - c0; acc = np.zeros(len(sel), np.float32); cnt = np.zeros(len(sel), np.float32)
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    v = a[rr + dr, cc + dc]; m = np.isfinite(v); acc[m] += v[m]; cnt[m] += 1
            val[sel] = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
G['chm'] = val
# landslide history relative to the shot time
CY = pd.read_parquet('v1/cell_year_events.parquet'); CY = CY[CY.cell.isin(G.cell[G.a > 0].unique())]
m = G[G.a > 0][['shot_number', 'cell', 't']].merge(CY[['cell', 't_event', 't_after', 't_before', 'frac']], on='cell')
m['before'] = m.t_after <= m.t; m['after'] = m.t_before > m.t; m['amb'] = ~(m.before | m.after)
h = m.groupby('shot_number').agg(n_before=('before', 'sum'), n_after=('after', 'sum'), n_amb=('amb', 'sum'),
                                 t_last_before=('t_event', lambda v: np.nan), ).reset_index()
tb = m[m.before].groupby('shot_number').agg(t_last=('t_event', 'max'), frac_before=('frac', 'max'))
ta = m[m.after].groupby('shot_number').agg(t_next=('t_event', 'min'), frac_after=('frac', 'max'))
h = h.drop(columns='t_last_before').merge(tb, on='shot_number', how='left').merge(ta, on='shot_number', how='left')
G = G.merge(h, on='shot_number', how='left')
G.to_parquet('v2/gedi_l4a_joined.parquet')
print(len(G), 'shots on grid;', int(G.chm.notna().sum()), 'with CHM;', int((G.a > 0).sum()), 'on catalogued cells;', int(((G.a >= 0.8) & (G.a3 >= 0.5)).sum()), 'on scar interiors')
print('typed forest shots', int(G.type.between(1, 6).sum()), 'never-catalogued within 5x5', int((G.type.between(1, 6) & ~G.ls5).sum()))
