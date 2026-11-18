# -*- coding: utf-8 -*-
# Rasterise the stand-volume map of the fourth national forest inventory (林地蓄積（碳儲存量）圖, MOA open data J18) to the same 30 m EPSG:3826 grid as forest_type_30m
# Usage: python rasterize_nfi4_volume.py <shapefile path without extension> <out .npz>
# (x0 151050, y1 2798910, 6680 x 12558). Layers: species (MajorTreeI code index), func (FunctionTy string index), volume (m3/ha, float32), area_ha.
import struct, numpy as np, cv2, json, sys, time
src = sys.argv[1]; out = sys.argv[2]; RES = 30.0; x0 = 151050.0; y1 = 2798910.0; W = 6680; H = 12558
SP = {'A134':'玉山圓柏','A114':'冷杉','A115':'雲杉','A101':'檜木（紅檜、扁柏）','A119':'台灣杉','A117':'杉木（巒大杉、杉木等）','A103':'臺灣肖楠','A118':'柳杉','A106':'松樹（二葉松、五葉松等）','A127':'南洋杉類','B620':'高山杜鵑類','B502':'相思樹','B301':'樟樹','B206':'大葉桃花心木','B515':'柚木','B507':'台灣赤楊','B504':'印度紫檀','B512':'桉樹類','B524':'銀合歡','B508':'木麻黃','B201':'烏心石','B516':'光臘樹','B202':'台灣櫸','B527':'木油桐','B630':'紅樹林類','C800':'單桿狀竹','A113':'鐵杉','C700':'叢生狀竹','B600':'其他(闊)','A151':'其他(針)'}
spidx = {k:i+1 for i,k in enumerate(SP)}
FN = ['1100','1200','1300','2100','2211','2212','2221','2222']; fnidx = {k:i+1 for i,k in enumerate(FN)}
f = open(src+'.dbf','rb'); hdr=f.read(32); nrec,hlen,rlen = struct.unpack('<IHH',hdr[4:12])
fields=[]
while True:
    fd=f.read(32)
    if fd[0]==0x0d: break
    fields.append((fd[:11].split(b'\0')[0].decode('latin1'),chr(fd[11]),fd[16]))
f.seek(hlen); sp=np.zeros(nrec,np.uint8); fn=np.zeros(nrec,np.uint8); vol=np.full(nrec,np.nan,np.float32)
for i in range(nrec):
    r=f.read(rlen); o=1; d={}
    for n,t,l in fields: d[n]=r[o:o+l]; o+=l
    sp[i]=spidx.get(d['MajorTreeI'].decode('utf-8','replace').strip(),0); fn[i]=fnidx.get(d['FunctionTy'].decode('latin1').strip(),0)
    try: vol[i]=float(d['Volume'])
    except: pass
s = open(src+'.shp','rb'); s.read(100)
S = np.zeros((H,W),np.uint8); Fm = np.zeros((H,W),np.uint8); V = np.full((H,W),np.nan,np.float32); t0=time.time(); i=0
while True:
    rh = s.read(8)
    if len(rh)<8: break
    recno, clen = struct.unpack('>ii', rh); body = s.read(clen*2); st = struct.unpack('<i', body[:4])[0]
    if st in (5,15,25):
        nparts, npts = struct.unpack('<ii', body[36:44]); parts = struct.unpack(f'<{nparts}i', body[44:44+4*nparts])
        pts = np.frombuffer(body[44+4*nparts:44+4*nparts+16*npts], dtype='<f8').reshape(npts,2)
        px = (pts[:,0]-x0)/RES; py = (y1-pts[:,1])/RES; bounds=list(parts)+[npts]
        polys=[np.round(np.stack([px[bounds[k]:bounds[k+1]], py[bounds[k]:bounds[k+1]]],1)).astype(np.int32) for k in range(nparts)]
        polys=[p for p in polys if len(p)>=3]
        if not polys: i+=1; continue
        if sp[i]: cv2.fillPoly(S, polys, int(sp[i]))
        if fn[i]: cv2.fillPoly(Fm, polys, int(fn[i]))
        if np.isfinite(vol[i]):
            # fill volume via a temporary uint8 mask restricted to the polygon bbox
            xs=np.concatenate([p[:,0] for p in polys]); ys=np.concatenate([p[:,1] for p in polys])
            xa,xb=max(xs.min(),0),min(xs.max()+1,W); ya,yb=max(ys.min(),0),min(ys.max()+1,H)
            if xb>xa and yb>ya:
                m=np.zeros((yb-ya,xb-xa),np.uint8); cv2.fillPoly(m,[p-np.array([xa,ya]) for p in polys],1); sub=V[ya:yb,xa:xb]; sub[m==1]=vol[i]
    i+=1
    if i%50000==0: print(i, f'{time.time()-t0:.0f}s', flush=True)
print('records', i, flush=True)
np.savez_compressed(out, species=S, func=Fm, volume=V)
json.dump(dict(crs='EPSG:3826', res=RES, x0=x0, y1=y1, width=W, height=H, transform=[RES,0,x0,0,-RES,y1], species_codes={v:{'code':k,'name':SP[k]} for k,v in spidx.items()}, func_codes={v:k for k,v in fnidx.items()},
               func_meaning={'1100':'原生林','1200':'原生林(保護)','1300':'原生林','2100':'半天然林','2211':'生產性人工林','2212':'保護性人工林','2221':'經改造天然林','2222':'經改造天然林'}, volume_unit='m3/ha (Volume field, 公頃材積)', n_records=int(i),
               source='Dataset/01_SOURCE/Forestry_Reference/Taiwan_4th_inventory_volume_1061/第四次森林資源調查林地蓄積圖1061.shp', note='same grid as Taiwan_forest_type_map_30m_TWD97/forest_type_30m; func_meaning for 1200/1300/2222 inferred from FunctionTy.xls not readable here — verify'), open(out+'.json','w'), ensure_ascii=False, indent=1)
u,c=np.unique(S,return_counts=True); print({int(k):round(float(v)*0.09) for k,v in zip(u,c)}, 'ha per species idx'); print('volume valid ha', round(float(np.isfinite(V).sum())*0.09), 'median', float(np.nanmedian(V)))
