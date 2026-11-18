# -*- coding: utf-8 -*-
"""s05: extend the cell table (stand origin, Global Forest Change loss year, polygon size, interior index); writes v2/cells.parquet."""
import os, sys, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
from config import SRC as DS, CUR, WORK
Q1 = WORK + '/q01'
X0, Y1, RES, W, H = 151050.0, 2798910.0, 30.0, 6680, 12558
A = pd.read_parquet('v1/cells.parquet'); CY = pd.read_parquet('v1/cell_year_events.parquet')
z = np.load(CUR + '/Forestry_Reference/Taiwan_forest_type_map_30m_TWD97/forest_type_30m.npz'); print(z.files)
FUNC = z['func']; T = z['type']
zv = np.load(CUR + '/Forestry_Reference/Taiwan_4th_NFI_volume_30m_TWD97/nfi4_volume_30m.npz'); print(zv.files)
def aligned(path, band):
    with rasterio.open(path) as ds:
        c0 = int(round((X0 - ds.bounds.left) / RES)); r0 = int(round((ds.bounds.top - Y1) / RES))
        a = ds.read(band, window=Window(c0, r0, W, H), boundless=True, fill_value=0)
    return a
with rasterio.open(DS + '/Earth_Observation/Lab_platform_assets/hansen_gfc_2025_30m.tif') as ds: print(ds.count, ds.descriptions, ds.dtypes)
LY = aligned(DS + '/Earth_Observation/Lab_platform_assets/hansen_gfc_2025_30m.tif', 2); GN = aligned(DS + '/Earth_Observation/Lab_platform_assets/hansen_gfc_2025_30m.tif', 3)
c = A.cell.values
A['func'] = FUNC.ravel()[c]; A['func_v'] = zv['func'].ravel()[c]; A['species'] = zv['species'].ravel()[c]
A['lossyear'] = np.nan_to_num(LY.ravel()[c]).astype(np.int16); A['gfc_gain'] = np.nan_to_num(GN.ravel()[c]).astype(np.int8)
# landslide-fraction raster and its 3x3 mean (interior index); never-catalogued neighbourhood for reference cells
LS = np.zeros(H * W, np.float32); m = ~A.is_ref.values; LS[c[m]] = A.a.values[m]; LS = LS.reshape(H, W)
P = np.pad(LS, 1); S = sum(P[1 + dr:1 + dr + H, 1 + dc:1 + dc + W] for dr in (-1, 0, 1) for dc in (-1, 0, 1)) / 9.0
A['a3'] = S.ravel()[c].astype(np.float32)
# largest polygon touching the cell (ha)
zz = np.load(Q1 + '/q01_record_cells.npz'); rec = pd.read_csv(Q1 + '/q01_records.csv', usecols=['fid', 'area_best_ha'])
pa = pd.DataFrame({'fid': zz['fid'], 'cell': zz['cell'].astype(np.int64)}).merge(rec, on='fid').groupby('cell').area_best_ha.max()
A['poly_ha'] = pa.reindex(c).values.astype(np.float32)
A.to_parquet('v2/cells.parquet'); np.save('v2/ls_frac_30m.npy', LS)
print(A[~A.is_ref].func.value_counts().sort_index().to_dict(), '| ref', A[A.is_ref].func.value_counts().sort_index().to_dict())
L = A[~A.is_ref]
print('lossyear>0 share (landslide cells, tc2000>=50):', round(float((L[L.tc2000 >= 50].lossyear > 0).mean()), 3), '| ref:', round(float((A[A.is_ref & (A.tc2000 >= 50)].lossyear > 0).mean()), 4))
core = L[(L.a >= 0.8) & (L.nyears == 1) & (L.tc2000 >= 50)]
d = (2000 + core.lossyear.where(core.lossyear > 0)) - core['first']
print('core single-event cells', len(core), 'flagged by GFC', round(float(core.lossyear.gt(0).mean()), 3), 'within -1..+1 yr of catalogue year', round(float(d.between(-1, 1).mean()), 3), 'within 0..+2', round(float(d.between(0, 2).mean()), 3))
print(d.value_counts().sort_index().loc[-5:5].to_dict())
