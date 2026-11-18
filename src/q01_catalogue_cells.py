# -*- coding: utf-8 -*-
"""q01: event-based landslide catalogue 2004-2025 (curated GeoPackage) on the 30 m grid of the fourth national forest inventory rasters.
Pure Python (sqlite3 + struct + cv2); GDAL is not required.
Each polygon is rasterised on a 6 m sub-grid (5x5 sub-cells of the 30 m NFI4 grid) so that small scars are
area-weighted instead of being snapped to whole 30 m cells.
Outputs: per-record table, per-(record, 30 m cell) fractions, yearly summaries.
"""
import sqlite3, struct, json, os, sys, time
import numpy as np, cv2, pandas as pd
from matplotlib.path import Path

from config import CUR
GPKG = CUR + '/Hazard_Events/Landslide_Taiwan/Phase19_Event_Catalog_2004_2025_20260825/geodata/phase19_event_catalog_candidates.gpkg'
FT = CUR + '/Forestry_Reference/Taiwan_forest_type_map_30m_TWD97/forest_type_30m.npz'
NV = CUR + '/Forestry_Reference/Taiwan_4th_NFI_volume_30m_TWD97/nfi4_volume_30m.npz'
OUT = sys.argv[1]
X0, Y1, RES, W, H = 151050.0, 2798910.0, 30.0, 6680, 12558
SUB = 5; RF = RES / SUB; SHIFT = 4; SC = 1 << SHIFT

def parse_gpkg(b):
    flags = b[3]; env = (flags >> 1) & 7
    off = 8 + {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]
    polys = []
    def read_geom(o):
        bo = '<' if b[o] == 1 else '>'
        gt = struct.unpack(bo + 'I', b[o+1:o+5])[0]; o += 5
        dim = 2 + (1 if (gt // 1000) in (1, 3) else 0) + (1 if (gt // 1000) in (2, 3) else 0)
        base = gt % 1000
        if base == 3:
            nr = struct.unpack(bo + 'I', b[o:o+4])[0]; o += 4; rings = []
            for _ in range(nr):
                npt = struct.unpack(bo + 'I', b[o:o+4])[0]; o += 4
                a = np.frombuffer(b, dtype=bo + 'f8', count=npt * dim, offset=o).reshape(npt, dim)[:, :2]
                o += npt * dim * 8; rings.append(a)
            polys.append(rings)
        elif base in (6, 7):
            ng = struct.unpack(bo + 'I', b[o:o+4])[0]; o += 4
            for _ in range(ng): o = read_geom(o)
        else:
            raise ValueError('geom type %d' % gt)
        return o
    read_geom(off)
    return polys

t0 = time.time()
T = np.load(FT)['type']; z = np.load(NV); V = z['volume']; SP = z['species']; VF = z['func']
print('rasters loaded', T.shape, V.shape, round(time.time() - t0, 1), flush=True)
con = sqlite3.connect('file:' + GPKG + '?mode=ro', uri=True)
cur = con.execute('select fid, record_id, catalog_year, event_category, area_best_ha, before_date, after_date, geom from event_catalogs')
rec_rows = []; cell_fid = []; cell_idx = []; cell_frac = []
nbad = 0
for k, (fid, rid, yr, cat, aha, bd, ad, g) in enumerate(cur):
    try:
        polys = parse_gpkg(g)
    except Exception as e:
        nbad += 1; continue
    allpts = np.vstack([r for p in polys for r in p])
    fx0 = int(np.floor((allpts[:, 0].min() - X0) / RF)) - 1; fx1 = int(np.ceil((allpts[:, 0].max() - X0) / RF)) + 1
    fy0 = int(np.floor((Y1 - allpts[:, 1].max()) / RF)) - 1; fy1 = int(np.ceil((Y1 - allpts[:, 1].min()) / RF)) + 1
    w = fx1 - fx0 + 1; h = fy1 - fy0 + 1
    m = np.zeros((h, w), np.uint8)
    for rings in polys:
        # 1) cv2 fill gives a superset (boundary-touching sub-cells included); 2) exact centre-in-polygon test (even-odd, holes respected)
        mm = np.zeros((h, w), np.uint8)
        pts = [np.round((np.stack([(r[:, 0] - X0) / RF - 0.5 - fx0, (Y1 - r[:, 1]) / RF - 0.5 - fy0], 1)) * SC).astype(np.int32) for r in rings]
        cv2.fillPoly(mm, pts, 1, lineType=cv2.LINE_8, shift=SHIFT)
        yy, xx = np.nonzero(mm)
        if len(yy) == 0: continue
        verts = []; codes = []
        for r in rings:
            verts.append(r); codes += [Path.MOVETO] + [Path.LINETO] * (len(r) - 2) + [Path.CLOSEPOLY]
        path = Path(np.vstack(verts), codes)
        cxy = np.stack([X0 + (xx + fx0 + 0.5) * RF, Y1 - (yy + fy0 + 0.5) * RF], 1)
        ins = path.contains_points(cxy)
        m[yy[ins], xx[ins]] = 1
    ys, xs = np.nonzero(m)
    if len(ys) == 0:
        rec_rows.append((fid, rid, yr, cat, aha, bd, ad, 0.0)); continue
    cy = (ys + fy0) // SUB; cx = (xs + fx0) // SUB
    ok = (cy >= 0) & (cy < H) & (cx >= 0) & (cx < W)
    flat = cy[ok].astype(np.int64) * W + cx[ok]
    u, c = np.unique(flat, return_counts=True)
    cell_fid.append(np.full(len(u), fid, np.int32)); cell_idx.append(u); cell_frac.append((c / (SUB * SUB)).astype(np.float32))
    rec_rows.append((fid, rid, yr, cat, aha, bd, ad, len(ys) * RF * RF / 1e4))
    if k % 5000 == 0: print(k, round(time.time() - t0, 1), flush=True)
print('rasterised', len(rec_rows), 'bad', nbad, round(time.time() - t0, 1), flush=True)
rec = pd.DataFrame(rec_rows, columns=['fid', 'record_id', 'catalog_year', 'event_category', 'area_best_ha', 'before_date', 'after_date', 'area_raster_ha'])
cf = np.concatenate(cell_fid); ci = np.concatenate(cell_idx); fr = np.concatenate(cell_frac)
cells = pd.DataFrame({'fid': cf, 'cell': ci, 'frac': fr})
cells['type'] = T.ravel()[ci]; cells['vol'] = V.ravel()[ci]; cells['species'] = SP.ravel()[ci]; cells['vfunc'] = VF.ravel()[ci]
cells = cells.merge(rec[['fid', 'catalog_year', 'event_category']], on='fid')
cells['ha'] = cells['frac'] * RES * RES / 1e4
os.makedirs(OUT, exist_ok=True)
rec.to_csv(OUT + '/q01_records.csv', index=False, encoding='utf-8-sig')
np.savez_compressed(OUT + '/q01_record_cells.npz', fid=cells['fid'].values, cell=cells['cell'].values, frac=cells['frac'].values,
                    year=cells['catalog_year'].values.astype(np.int16), type=cells['type'].values, vol=cells['vol'].values,
                    species=cells['species'].values, vfunc=cells['vfunc'].values)
# national context of the two NFI4 rasters
ctx = {'type_ha': {int(k): float(v) * RES * RES / 1e4 for k, v in zip(*np.unique(T, return_counts=True))},
       'volume_cells_gt0_ha': float((V > 0).sum()) * 0.09, 'volume_total_m3': float(V[V > 0].sum() * 0.09),
       'mapped_volume_layer_ha': float((SP > 0).sum()) * 0.09}
for t in range(1, 10):
    mk = (T == t); vv = V[mk]
    ctx['type%d' % t] = {'ha': float(mk.sum()) * 0.09, 'ha_vol_gt0': float((vv > 0).sum()) * 0.09, 'mean_vol_gt0': float(vv[vv > 0].mean()) if (vv > 0).any() else None}
json.dump(ctx, open(OUT + '/q01_nfi4_context.json', 'w'), ensure_ascii=False, indent=1)
# yearly summaries
g = cells.groupby(['catalog_year', 'type'])['ha'].sum().unstack(fill_value=0.0)
g.to_csv(OUT + '/q01_year_by_type_ha.csv', encoding='utf-8-sig')
cells['volha'] = cells['vol'] * cells['ha']
cells['has_vol'] = (cells['vol'] > 0)
s = cells.groupby('catalog_year').apply(lambda d: pd.Series({'ha_total': d['ha'].sum(), 'ha_vol_gt0': d.loc[d.has_vol, 'ha'].sum(), 'vol_m3': d['volha'].sum()}))
s.to_csv(OUT + '/q01_year_volume.csv', encoding='utf-8-sig')
gv = cells[cells.has_vol].groupby(['catalog_year', 'type']).apply(lambda d: pd.Series({'ha': d['ha'].sum(), 'vol_m3': d['volha'].sum()}))
gv.to_csv(OUT + '/q01_year_type_volume.csv', encoding='utf-8-sig')
chk = rec.groupby('catalog_year')[['area_best_ha', 'area_raster_ha']].sum()
chk.to_csv(OUT + '/q01_area_check.csv', encoding='utf-8-sig')
print(chk.round(1).to_string()); print('DONE', round(time.time() - t0, 1), flush=True)
