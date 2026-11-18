# -*- coding: utf-8 -*-
# Rasterisation of the forest-type map of the fourth national forest inventory (全島森林林型分布圖, data.gov.tw dataset 57873) onto the
# project 30 m grid in TWD97 / TM2 (EPSG:3826). Pure Python + OpenCV; GDAL is not required.
# Usage: python rasterize_forest_type.py <shapefile path without extension> <out .npz>
import struct, numpy as np, cv2, json, sys, time
src = sys.argv[1]; out = sys.argv[2]; RES = 30.0
TYPES = {'闊葉樹林型':1,'針葉樹林型':2,'竹林':3,'針闊葉樹混淆':4,'針闊葉樹混':4,'竹闊混淆林':5,'竹針混淆林':6,'竹針闊混淆林':6,'待成林地':7,'裸露地':8,'陰影':9}
FUNC = {'原生林':1,'經改造天然林':2,'半天然林':3,'生產性人工林':4,'保護性人工林':5,'':0}
# dbf
f = open(src+'.dbf','rb'); hdr=f.read(32); nrec,hlen,rlen = struct.unpack('<IHH',hdr[4:12]); f.seek(hlen)
tcode = np.zeros(nrec, np.uint8); fcode = np.zeros(nrec, np.uint8)
for i in range(nrec):
    rec=f.read(rlen); tcode[i]=TYPES.get(rec[1:21].decode('utf-8','replace').strip(),0); fcode[i]=FUNC.get(rec[31:81].decode('utf-8','replace').strip(),0)
# shp header bbox
s = open(src+'.shp','rb'); h = s.read(100)
xmin,ymin,xmax,ymax = struct.unpack('<4d', h[36:68]); print('bbox', xmin,ymin,xmax,ymax, flush=True)
x0 = np.floor(xmin/RES)*RES; y1 = np.ceil(ymax/RES)*RES
W = int(np.ceil((xmax-x0)/RES)); H = int(np.ceil((y1-ymin)/RES)); print('grid', W, H, flush=True)
T = np.zeros((H,W), np.uint8); F = np.zeros((H,W), np.uint8)
t0=time.time(); i=0
while True:
    rh = s.read(8)
    if len(rh)<8: break
    recno, clen = struct.unpack('>ii', rh); body = s.read(clen*2)
    st = struct.unpack('<i', body[:4])[0]
    if st in (5,15,25):
        nparts, npts = struct.unpack('<ii', body[36:44])
        parts = struct.unpack(f'<{nparts}i', body[44:44+4*nparts])
        pts = np.frombuffer(body[44+4*nparts:44+4*nparts+16*npts], dtype='<f8').reshape(npts,2)
        px = ((pts[:,0]-x0)/RES); py = ((y1-pts[:,1])/RES)
        polys=[]; bounds=list(parts)+[npts]
        for k in range(nparts):
            seg = np.stack([px[bounds[k]:bounds[k+1]], py[bounds[k]:bounds[k+1]]],1)
            polys.append(np.round(seg).astype(np.int32))
        tc = int(tcode[i]); fc = int(fcode[i])
        if tc: cv2.fillPoly(T, polys, tc)   # holes: shapefile rings are CW outer/CCW inner; fillPoly with all rings uses even-odd -> holes handled
        if fc: cv2.fillPoly(F, polys, fc)
    i += 1
    if i % 50000 == 0: print(i, f'{time.time()-t0:.0f}s', flush=True)
print('records', i, flush=True)
np.savez_compressed(out, type=T, func=F)
json.dump(dict(crs='EPSG:3826', res=RES, x0=float(x0), y1=float(y1), width=W, height=H, transform=[RES,0,float(x0),0,-RES,float(y1)],
               type_codes={v:k for k,v in TYPES.items() if v!=4 or k=='針闊葉樹混淆'}, func_codes={v:k for k,v in FUNC.items()}, n_records=int(i),
               source='Dataset/01_SOURCE/Forestry_Reference/Taiwan_forest_type_map/全島森林林型分布圖.shp (4th NFI)'), open(out+'.json','w'), ensure_ascii=False, indent=1)
u,c = np.unique(T, return_counts=True); print(dict(zip(u.tolist(), (c*RES*RES/1e4).round().tolist())), 'ha per class')
