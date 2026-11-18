"""NDVI time series of the Laonong River case site (Fig. 7e): Landsat Collection 2 Level-2 scenes of years y0..y1 with cloud cover < 35%,
at most four clearest scenes per year with more than half of the window clear; per scene the median NDVI of each cell group.
Usage (from the working directory): python src/cases/ts.py 2005 2012"""
import sys,time,json,os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from stac import *; from chip import chip
import numpy as np, pandas as pd
R=(7500,7800,2300,2700); r0,r1,c0,c1=R; W=6680
y0,y1=int(sys.argv[1]),int(sys.argv[2])
C=pd.read_parquet('v2/cells.parquet'); C['row']=C.cell//W; C['col']=C.cell%W
c=C[(C.row>=r0)&(C.row<r1)&(C.col>=c0)&(C.col<c1)]
G={}
s=c[(c.is_ref==0)&(c['first']==2009)&(c.a>=0.8)]
G['large']=s[s.poly_ha>=100]; G['small']=s[(s.poly_ha<3)&(s.nyears==1)&(s.tc2000>=75)]
G['ref']=c[c.is_ref==1]
idx={k:((v.row-r0).values,(v.col-c0).values) for k,v in G.items()}
bb=bbox_ll(*R); out=[]
for yr in range(y0,y1+1):
    L=search('landsat-c2-l2',bb,f'{yr}-01-01/{yr}-12-31',{'eo:cloud_cover':{'lt':35}})
    L=[f for f in L if f['properties']['platform']!='landsat-7' or yr<2003 or True]
    L=sorted(L,key=lambda f:f['properties']['eo:cloud_cover'])
    n=0
    for f in L:
        if n>=4: break
        try: b,m=chip(f,*R,30)
        except Exception as ex: print('ERR',f['id'],ex); continue
        if m.mean()<0.5: continue
        nd=(b[3]-b[2])/(b[3]+b[2]); nd[~m]=np.nan
        rec={'id':f['id'],'date':f['properties']['datetime'][:10],'platform':f['properties']['platform'],'clear':float(m.mean())}
        for k,(rr,cc) in idx.items():
            v=nd[rr,cc]; rec[k]=float(np.nanmedian(v)) if np.isfinite(v).sum()>=0.5*len(v) else np.nan; rec[k+'_n']=int(np.isfinite(v).sum())
        out.append(rec); n+=1; print(rec['date'],rec['platform'],round(rec['clear'],2),{k:round(rec[k],3) for k in G},flush=True)
os.makedirs('cases', exist_ok=True); pd.DataFrame(out).to_csv(f'cases/A_ndvi_{y0}_{y1}.csv',index=False)
print({k:len(v) for k,v in G.items()})
