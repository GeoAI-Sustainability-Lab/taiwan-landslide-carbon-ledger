# -*- coding: utf-8 -*-
"""CHMv2 observation-date polygons (53 tile GeoJSONs, WGS84) -> 30 m TWD97 raster of acquisition date
(uint16 days since 2000-01-01; 0 = no date). Pure numpy + cv2; GDAL and pyproj are not required."""
import json, glob, os, sys, time, datetime
import numpy as np, cv2
from config import SRC as _SRC
SRC = _SRC + '/Terrain_Canopy/CHMv2_Meta_WRI_2026/metadata/tile_observation_dates'
OUT = sys.argv[1]
X0, Y1, RES, W, H = 151050.0, 2798910.0, 30.0, 6680, 12558
a = 6378137.0; f = 1 / 298.257222101; k0 = 0.9999; lon0 = np.radians(121.0); FE = 250000.0
e2 = f * (2 - f); ep2 = e2 / (1 - e2)
def tm(lon, lat):
    ph = np.radians(lat); la = np.radians(lon)
    N = a / np.sqrt(1 - e2 * np.sin(ph) ** 2); T = np.tan(ph) ** 2; C = ep2 * np.cos(ph) ** 2; A = (la - lon0) * np.cos(ph)
    M = a * ((1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256) * ph - (3 * e2 / 8 + 3 * e2 ** 2 / 32 + 45 * e2 ** 3 / 1024) * np.sin(2 * ph)
             + (15 * e2 ** 2 / 256 + 45 * e2 ** 3 / 1024) * np.sin(4 * ph) - (35 * e2 ** 3 / 3072) * np.sin(6 * ph))
    x = FE + k0 * N * (A + (1 - T + C) * A ** 3 / 6 + (5 - 18 * T + T ** 2 + 72 * C - 58 * ep2) * A ** 5 / 120)
    y = k0 * (M + N * np.tan(ph) * (A ** 2 / 2 + (5 - T + 9 * C + 4 * C ** 2) * A ** 4 / 24 + (61 - 58 * T + T ** 2 + 600 * C - 330 * ep2) * A ** 6 / 720))
    return x, y
R = np.zeros((H, W), np.uint16); t0 = time.time(); base = datetime.date(2000, 1, 1); nfeat = 0
for p in sorted(glob.glob(SRC + '/*.geojson')):
    fc = json.load(open(p))
    for ft in fc['features']:
        d = ft['properties'].get('acq_date')
        if not d: continue
        v = (datetime.date.fromisoformat(d[:10]) - base).days
        polys = []
        def collect(g):
            if g is None: return
            if g['type'] == 'Polygon': polys.append(g['coordinates'])
            elif g['type'] == 'MultiPolygon': polys.extend(g['coordinates'])
            elif g['type'] == 'GeometryCollection':
                for q in g['geometries']: collect(q)
        collect(ft['geometry'])
        for rings in polys:
            pts = []
            for r in rings:
                r = np.asarray(r, float); x, y = tm(r[:, 0], r[:, 1])
                pts.append(np.round(np.stack([(x - X0) / RES - 0.5, (Y1 - y) / RES - 0.5], 1)).astype(np.int32))
            cv2.fillPoly(R, pts, int(v))
        nfeat += 1
    del fc
    print(os.path.basename(p), nfeat, round(time.time() - t0), flush=True)
np.savez_compressed(OUT, days=R)
yrs = 2000 + R[R > 0] / 365.25
print('DONE cells with date', int((R > 0).sum()), 'year quantiles', np.round(np.quantile(yrs, [0.05, 0.25, 0.5, 0.75, 0.95]), 2), flush=True)
