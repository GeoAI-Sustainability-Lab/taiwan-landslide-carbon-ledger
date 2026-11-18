# -*- coding: utf-8 -*-
"""Compile the downloaded MapServer pages of layer 1 (third national forest inventory photo points; download_nfi3_photo_points.py)
into the table read by s14_nfi3_check.py: OBJECTID, geom_x, geom_y (TWD97, as served), DATE, LANDUSE, AGE, DENSITY, STAND, VOLUME,
ELEVATION, OWNERSHIP. Usage: python compile_nfi3_photo_points.py <page_dir> <out_csv>"""
import sys, glob, json
import pandas as pd
KEEP = ['DATE', 'LANDUSE', 'AGE', 'DENSITY', 'STAND', 'VOLUME', 'ELEVATION', 'OWNERSHIP']
rows = []
for fn in sorted(glob.glob(sys.argv[1] + '/page_*.json')):
    for f in json.load(open(fn, encoding='utf8'))['features']:
        a = f['attributes']; g = f.get('geometry') or {}
        rows.append({'OBJECTID': a['OBJECTID'], 'geom_x': g.get('x'), 'geom_y': g.get('y'), **{k: a.get(k) for k in KEEP}})
pd.DataFrame(rows).sort_values('OBJECTID').to_csv(sys.argv[2], index=False)
print(len(rows), 'points ->', sys.argv[2])
