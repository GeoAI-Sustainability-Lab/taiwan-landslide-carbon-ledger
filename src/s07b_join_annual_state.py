# -*- coding: utf-8 -*-
"""s07b: add the official annual landslide maps 2017-2024 (s07 output annual_state/annual_state_<year>.npz: cell, sub, year) to the
cell table as the fraction of each 30 m cell mapped as landslide (columns s2017 ... s2024), and the forest-type group tg.
Run from the working directory after s05."""
import glob, numpy as np, pandas as pd
from config import W, H
TG = {1: 'BL', 2: 'CF', 4: 'MX', 3: 'BB', 5: 'BB', 6: 'BB'}
A = pd.read_parquet('v2/cells.parquet'); A['tg'] = A.type_nb.map(TG)
for f in sorted(glob.glob('annual_state/annual_state_*.npz')):
    z = np.load(f); y = int(z['year']); r = np.zeros(W * H, np.uint8); r[z['cell']] = z['sub']
    A[f's{y}'] = r[A.cell.values] / 72.0
A.to_parquet('v2/cells.parquet')
print('annual maps joined:', sorted(c for c in A.columns if c.startswith('s20')))
