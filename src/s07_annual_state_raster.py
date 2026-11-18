# -*- coding: utf-8 -*-
"""Rasterise an official annual whole-island landslide layer (ARDSWC, polygons >= 0.1 ha) onto the project 30 m TWD97 grid.
Each 30 m cell stores the number of half 5 m sub-cells (0-72; boundary sub-cells count one half) inside mapped landslide polygons. Pure numpy + OpenCV."""
import sys, os, zipfile, struct, time, numpy as np, cv2
y = int(sys.argv[1]); src = sys.argv[2]; out = sys.argv[3]
X0, Y1, RES, W, H, SUB = 151050.0, 2798910.0, 30.0, 6680, 12558, 6
t0 = time.time()
if src.endswith('.zip'):
    z = zipfile.ZipFile(src); nm = [n for n in z.namelist() if n.lower().endswith('.shp')][0]; b = z.read(nm)
else: b = open(src, 'rb').read()
print('read', len(b), round(time.time() - t0, 1), flush=True)
cnt = np.zeros((H, W), np.uint8); pos = 100; n = 0; nbig = 0; area = 0.0; d = RES / SUB; SH = 4; F = 1 << SH
while pos < len(b):
    clen = struct.unpack('>i', b[pos + 4:pos + 8])[0]; p = pos + 8; pos = p + clen * 2
    st = struct.unpack('<i', b[p:p + 4])[0]
    if st not in (5, 15, 25): continue
    xmin, ymin, xmax, ymax = struct.unpack('<4d', b[p + 4:p + 36]); nparts, npts = struct.unpack('<ii', b[p + 36:p + 44])
    parts = np.frombuffer(b, '<i4', nparts, p + 44); pts = np.frombuffer(b, '<f8', npts * 2, p + 44 + 4 * nparts).reshape(-1, 2)
    c0 = max(int(np.floor((xmin - X0) / RES)), 0); c1 = min(int(np.ceil((xmax - X0) / RES)), W); r0 = max(int(np.floor((Y1 - ymax) / RES)), 0); r1 = min(int(np.ceil((Y1 - ymin) / RES)), H)
    if c1 <= c0 or r1 <= r0: continue
    n += 1
    # sub-pixel coordinates with pixel centres at integers, fixed-point with SH fractional bits
    px = (pts[:, 0] - (X0 + c0 * RES)) / d - 0.5; py = ((Y1 - r0 * RES) - pts[:, 1]) / d - 0.5
    idx = list(parts) + [npts]
    STR = 600   # 30 m rows per strip
    if (r1 - r0) > STR or (c1 - c0) > 3000: nbig += 1
    for s0 in range(r0, r1, STR):
        s1 = min(s0 + STR, r1); hh = (s1 - s0) * SUB; ww = (c1 - c0) * SUB
        m = np.zeros((hh, ww), np.uint8)
        rings = [np.round(np.column_stack([px[idx[k]:idx[k + 1]], py[idx[k]:idx[k + 1]] - (s0 - r0) * SUB]) * F).astype(np.int32) for k in range(nparts)]
        cv2.fillPoly(m, rings, 2, lineType=cv2.LINE_8, shift=SH); e = np.zeros((hh, ww), np.uint8); cv2.polylines(e, rings, True, 1, lineType=cv2.LINE_8, shift=SH)
        m = np.where((e > 0) & (m > 0), 1, m)   # boundary pixels count one half
        a = m.reshape(s1 - s0, SUB, c1 - c0, SUB).sum(axis=(1, 3))
        v = cnt[s0:s1, c0:c1].astype(np.int32) + a; cnt[s0:s1, c0:c1] = np.minimum(v, 72).astype(np.uint8)
    # shoelace area (m2), holes negative by ring orientation
    for k in range(nparts):
        q = pts[idx[k]:idx[k + 1]]; area += -0.5 * float(np.sum(q[:-1, 0] * q[1:, 1] - q[1:, 0] * q[:-1, 1]))
nz = np.flatnonzero(cnt)
np.savez_compressed(out, cell=nz.astype(np.int32), sub=cnt.ravel()[nz], year=y)
print(y, 'polygons', n, 'large', nbig, 'vector ha', round(abs(area) / 1e4, 1), 'raster ha', round(float(cnt.sum()) / 2 * d * d / 1e4, 1), 'cells', len(nz), 'sec', round(time.time() - t0, 1), flush=True)
