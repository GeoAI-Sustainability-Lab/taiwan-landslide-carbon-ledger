"""Image chips on the project grid from Planetary Computer STAC items: Landsat Collection 2 Level-2 and Sentinel-2 Level-2A surface reflectance
(blue, green, red, NIR, SWIR1) with a clear-sky mask from the QA_PIXEL or SCL band."""
import os, numpy as np, rasterio, sys
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from stac import sign
X0,Y1=151050.0,2798910.0
LB={'landsat-5':['blue','green','red','nir08','swir16'],'landsat-7':['blue','green','red','nir08','swir16'],
    'landsat-8':['blue','green','red','nir08','swir16'],'landsat-9':['blue','green','red','nir08','swir16']}
SB=['B02','B03','B04','B08','B11']
def _rd(href,coll,dst,dst_t,rs):
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_MAX_RETRY='4',GDAL_HTTP_RETRY_DELAY='2'):
        with rasterio.open(sign(href,coll)) as s:
            reproject(rasterio.band(s,1),dst,dst_transform=dst_t,dst_crs='EPSG:3826',resampling=rs,src_nodata=s.nodata,dst_nodata=0)
def chip(item,r0,r1,c0,c1,res):
    """returns dict: bands (5,h,w) reflectance [blue,green,red,nir,swir1], clear mask"""
    f=int(30/res); h,w=(r1-r0)*f,(c1-c0)*f; dst_t=from_origin(X0+c0*30,Y1-r0*30,res,res)
    coll=item['collection']; A=item['assets']
    if coll=='landsat-c2-l2':
        names=LB[item['properties']['platform']]; out=np.zeros((5,h,w),np.float32)
        for i,n in enumerate(names):
            d=np.zeros((h,w),np.float32); _rd(A[n]['href'],coll,d,dst_t,Resampling.bilinear); out[i]=np.where(d>0,d*2.75e-5-0.2,np.nan)
        q=np.zeros((h,w),np.float32); _rd(A['qa_pixel']['href'],coll,q,dst_t,Resampling.nearest); q=q.astype(np.uint32)
        bad=(q==0)|((q>>1)&1).astype(bool)|((q>>3)&1).astype(bool)|((q>>4)&1).astype(bool)|((q>>2)&1).astype(bool)
    else:
        out=np.zeros((5,h,w),np.float32)
        base=item['properties'].get('s2:processing_baseline','0'); off=-1000 if float(base)>=4.0 else 0
        for i,n in enumerate(SB):
            d=np.zeros((h,w),np.float32); _rd(A[n]['href'],coll,d,dst_t,Resampling.bilinear); out[i]=np.where(d>0,(d+off)/1e4,np.nan)
        q=np.zeros((h,w),np.float32); _rd(A['SCL']['href'],coll,q,dst_t,Resampling.nearest)
        bad=np.isin(q,[0,1,3,8,9,10])
    bad|=np.isnan(out).any(0)
    return out,~bad
