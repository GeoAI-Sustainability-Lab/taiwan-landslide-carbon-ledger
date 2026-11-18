"""CHMv2 10 m mean canopy height resampled onto a window of the project grid (rows r0:r1, columns c0:c1; 10 m pixels)."""
import os, sys, numpy as np, rasterio, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from config import SRC
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
W=6680; X0,Y1=151050.0,2798910.0
T = SRC + '/Terrain_Canopy/'
CHM=[T+'CHMv2_10m_mean_Z7_1321230_1321231_20261111/mean/1321230.tif',T+'CHMv2_10m_mean_Z7_1321230_1321231_20261111/mean/1321231.tif',
     T+'CHMv2_10m_mean_p95_Z7_1321232_1321233_20261107/mean/1321232.tif',T+'CHMv2_10m_mean_p95_Z7_1321232_1321233_20261107/mean/1321233.tif']
def chm10(r0,r1,c0,c1):
    dst_t=from_origin(X0+c0*30, Y1-r0*30, 10, 10); out=np.full(((r1-r0)*3,(c1-c0)*3),np.nan,np.float32)
    for p in CHM:
        d=np.full(out.shape,255,np.float32)
        with rasterio.open(p) as s:
            reproject(rasterio.band(s,1),d,dst_transform=dst_t,dst_crs='EPSG:3826',resampling=Resampling.nearest,src_nodata=255,dst_nodata=255)
        out=np.where(np.isnan(out)&(d!=255),d,out)
    return out
