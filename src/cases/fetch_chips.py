# -*- coding: utf-8 -*-
"""Image chips and canopy-height windows used in Figs 4 and 7 (run from the working directory; writes cases/).
Site A, Laonong River (rows 7500-7800, columns 2300-2700): Landsat 5 scenes of 24 Aug 2008 and 12 Sep 2009 (30 m) and Sentinel-2 scenes of
17 Nov 2017 and 5 Jan 2024 (10 m). Site B, Heping (rows 3650-3900, columns 5380-5640): Sentinel-2 scenes of 14 Feb 2024 and 25 Mar 2025 (10 m).
Each chip is saved as cases/<site>_<item id>.npz with b (blue, green, red, NIR, SWIR1 reflectance) and m (clear-sky mask)."""
import os, sys, numpy as np, requests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stac import API
from chip import chip
from common import chm10

RA, RB = (7500, 7800, 2300, 2700), (3650, 3900, 5380, 5640)
ITEMS = [('A', 'landsat-c2-l2', 'LT05_L2SP_117044_20080824_02_T1', 30), ('A', 'landsat-c2-l2', 'LT05_L2SP_117044_20090912_02_T1', 30),
         ('A', 'sentinel-2-l2a', 'S2B_MSIL2A_20171117T022939_R046_T50QRL_20201014T202927', 10), ('A', 'sentinel-2-l2a', 'S2B_MSIL2A_20240105T023109_R046_T50QRL_20240105T062216', 10),
         ('B', 'sentinel-2-l2a', 'S2B_MSIL2A_20240214T022819_R046_T51QUG_20240214T064747', 10), ('B', 'sentinel-2-l2a', 'S2C_MSIL2A_20250325T022551_R046_T51QUG_20250325T062212', 10)]

def item(coll, iid):
    r = requests.post(API, json={'collections': [coll], 'ids': [iid], 'limit': 1}, timeout=120); r.raise_for_status()
    f = r.json()['features']; assert f, iid
    return f[0]

if __name__ == '__main__':
    os.makedirs('cases', exist_ok=True)
    for site, coll, iid, res in ITEMS:
        out = f'cases/{site}_{iid}.npz'
        if os.path.exists(out): continue
        b, m = chip(item(coll, iid), *(RA if site == 'A' else RB), res)
        np.savez_compressed(out, b=b, m=m); print(out, b.shape, round(float(m.mean()), 3), flush=True)
    np.save('cases/A_chm10.npy', chm10(*RA)); np.save('cases/B_chm10.npy', chm10(*RB))
