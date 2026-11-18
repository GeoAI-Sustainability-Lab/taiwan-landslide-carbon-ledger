"""Planetary Computer STAC search and SAS-token signing."""
import requests, json, numpy as np, pyproj, time
API='https://planetarycomputer.microsoft.com/api/stac/v1/search'
_tok={}
def token(coll):
    if coll not in _tok:
        _tok[coll]=requests.get(f'https://planetarycomputer.microsoft.com/api/sas/v1/token/{coll}',timeout=60).json()['token']
    return _tok[coll]
def sign(href,coll): return href+('&' if '?' in href else '?')+token(coll)
TW=pyproj.Transformer.from_crs(3826,4326,always_xy=True)
def bbox_ll(r0,r1,c0,c1,X0=151050.0,Y1=2798910.0):
    xs=[X0+c0*30,X0+c1*30]; ys=[Y1-r1*30,Y1-r0*30]
    lo,la=TW.transform([xs[0],xs[1],xs[0],xs[1]],[ys[0],ys[0],ys[1],ys[1]])
    return [min(lo),min(la),max(lo),max(la)]
def search(coll,bbox,dt,query=None,limit=200):
    body={'collections':[coll],'bbox':bbox,'datetime':dt,'limit':limit}
    if query: body['query']=query
    out=[]; url=API
    while True:
        r=requests.post(url,json=body,timeout=120); r.raise_for_status(); j=r.json(); out+=j['features']
        nx=[l for l in j.get('links',[]) if l['rel']=='next']
        if not nx or len(out)>=1000: break
        body=nx[0].get('body',body); url=nx[0]['href']
    return out
