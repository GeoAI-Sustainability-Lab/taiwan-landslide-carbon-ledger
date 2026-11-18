"""Download a full ArcGIS MapServer layer page by page (raw JSON pages kept as provenance); layer 1 = third national forest inventory photo points.
Usage: python download_nfi3_photo_points.py <layer_id> <out_dir> [page_size]
Source: Forestry and Nature Conservation Agency public MapServer MA_Map0210.
"""
import sys, json, time, hashlib, os, datetime, urllib.parse
import requests

BASE = "https://gis.forest.gov.tw/arcgis/rest/services/MA_Map0210/MapServer"
layer = int(sys.argv[1]); out = sys.argv[2]
page = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
os.makedirs(out, exist_ok=True)
s = requests.Session()
s.headers["User-Agent"] = "NCHU-GeoAI-Sustainability-Lab research download (academic use)"

def get(url, params, tries=6):
    for k in range(tries):
        try:
            r = s.get(url, params=params, timeout=180)
            if r.status_code == 200 and r.content[:1] == b"{":
                j = r.json()
                if "error" not in j:
                    return r.content, j
                err = j["error"]
            else:
                err = f"http {r.status_code}"
        except Exception as e:
            err = repr(e)
        time.sleep(3 + 5 * k)
    raise RuntimeError(f"failed {params}: {err}")

meta_raw, meta = get(f"{BASE}/{layer}", {"f": "json"})
open(os.path.join(out, "layer_metadata.json"), "wb").write(meta_raw)
_, cnt = get(f"{BASE}/{layer}/query", {"where": "1=1", "returnCountOnly": "true", "f": "json"})
total = cnt["count"]
_, ids = get(f"{BASE}/{layer}/query", {"where": "1=1", "returnIdsOnly": "true", "f": "json"})
oids = sorted(ids["objectIds"])
oid_field = ids.get("objectIdFieldName", "OBJECTID")
print("layer", layer, meta.get("name"), "count", total, "ids", len(oids), flush=True)
receipt = {"service": BASE, "layer": layer, "name": meta.get("name"), "count_reported": total,
           "n_object_ids": len(oids), "page_size": page, "started_utc": datetime.datetime.utcnow().isoformat() + "Z",
           "pages": []}
n = 0
for i in range(0, len(oids), page):
    chunk = oids[i:i + page]
    fn = os.path.join(out, f"page_{i // page:04d}.json")
    if os.path.exists(fn) and os.path.getsize(fn) > 100:
        raw = open(fn, "rb").read(); j = json.loads(raw)
    else:
        where = f"{oid_field} >= {chunk[0]} AND {oid_field} <= {chunk[-1]}"
        raw, j = get(f"{BASE}/{layer}/query", {"where": where, "outFields": "*", "returnGeometry": "true",
                                               "orderByFields": oid_field, "f": "json"})
        open(fn, "wb").write(raw)
        time.sleep(0.4)
    k = len(j.get("features", []))
    n += k
    receipt["pages"].append({"file": os.path.basename(fn), "oid_min": chunk[0], "oid_max": chunk[-1], "features": k,
                             "expected": len(chunk), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                             "exceededTransferLimit": bool(j.get("exceededTransferLimit"))})
    if (i // page) % 10 == 0:
        print(f"page {i // page} cum {n}/{len(oids)}", flush=True)
receipt["features_total"] = n
receipt["finished_utc"] = datetime.datetime.utcnow().isoformat() + "Z"
receipt["complete"] = (n == len(oids))
json.dump(receipt, open(os.path.join(out, "download_receipt.json"), "w"), ensure_ascii=False, indent=1)
print("DONE", n, "of", len(oids), "complete", receipt["complete"], flush=True)
