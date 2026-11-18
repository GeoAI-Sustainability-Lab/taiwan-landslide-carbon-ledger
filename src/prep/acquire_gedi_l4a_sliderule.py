# -*- coding: utf-8 -*-
"""GEDI L4A footprints over Taiwan through the public SlideRule Earth service (no user credentials; the service subsets
ORNL DAAC granules server-side). Modes: cmr | test | run <out_dir> | decode <out_dir> <npz>"""
import sys, os, json, struct, time, hashlib, re
import requests, numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
URL = 'https://sliderule.slideruleearth.io/source/'
CMR = 'https://cmr.earthdata.nasa.gov/search/granules.json'
COLL = 'C2237824918-ORNL_CLOUD'      # GEDI L4A V2.1
BOX = (119.9, 21.8, 122.1, 25.4)
POLY = [{'lon': BOX[0], 'lat': BOX[1]}, {'lon': BOX[2], 'lat': BOX[1]}, {'lon': BOX[2], 'lat': BOX[3]}, {'lon': BOX[0], 'lat': BOX[3]}, {'lon': BOX[0], 'lat': BOX[1]}]
FMT = {'INT8': 'b', 'UINT8': 'B', 'INT16': 'h', 'UINT16': 'H', 'INT32': 'i', 'UINT32': 'I', 'INT64': 'q', 'UINT64': 'Q', 'TIME8': 'q', 'FLOAT': 'f', 'DOUBLE': 'd'}
def cmr():
    out = []; page = 1
    while True:
        r = requests.get(CMR, params={'collection_concept_id': COLL, 'bounding_box': '%g,%g,%g,%g' % BOX, 'page_size': 2000, 'page_num': page}, timeout=60); r.raise_for_status()
        e = r.json()['feed']['entry']
        if not e: break
        out += [{'id': x['id'], 'name': x.get('producer_granule_id') or x['title'], 'time_start': x.get('time_start')} for x in e]; page += 1
    return out
def defs():
    d = {}
    for k in ('gedi04arec', 'gedi04arec.footprint', 'exceptrec'):
        d[k] = requests.get(URL + 'definition', data=json.dumps({'rectype': k}), timeout=30).json()
    return d
def decode(data, D):
    foot = D['gedi04arec.footprint']; size = foot['__datasize']; off = D['gedi04arec']['footprint']['offset'] // 8
    dt = np.dtype({'names': [k for k in foot if not k.startswith('__')], 'formats': ['<' + FMT[foot[k]['type']] for k in foot if not k.startswith('__')],
                   'offsets': [foot[k]['offset'] // 8 for k in foot if not k.startswith('__')], 'itemsize': size})
    pos = 0; arrs = []; msgs = []
    while pos < len(data):
        ver, nt, nd = struct.unpack_from('>hhi', data, pos); end = pos + 8 + nt + nd; raw = data[pos + 8:end]
        kind = raw[:nt].rstrip(b'\0').decode('ascii'); payload = raw[nt:]
        if kind == 'gedi04arec': arrs.append(np.frombuffer(payload[off:off + ((len(payload) - off) // size) * size], dtype=dt).copy())
        elif kind == 'exceptrec':
            ex = D['exceptrec']; code = struct.unpack_from('<i', payload, ex['code']['offset'] // 8)[0]
            txt = payload[ex['text']['offset'] // 8:].split(b'\0')[0].decode('utf8', 'replace'); msgs.append((code, txt))
        pos = end
    return (np.concatenate(arrs) if arrs else np.zeros(0, dt)), msgs
def fetch(name, outdir=None):
    res = name.split('.', 1)[1] if name.startswith('GEDI_L4A') else name
    payload = {'resources': [res], 'parms': {'asset': 'gedil4a', 'poly': POLY, 'degrade_filter': True, 'l2_quality_filter': True, 'l4_quality_filter': True, 'surface_filter': True}}
    t0 = time.time(); r = requests.post(URL + 'gedi04ap', json=payload, timeout=(15, float(os.environ.get('RT', 600)))); body = r.content
    if outdir:
        open(os.path.join(outdir, res + '.bin'), 'wb').write(body)
    return res, r.status_code, len(body), round(time.time() - t0, 1), body
if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'cmr':
        g = cmr(); json.dump(g, open('granules_l4a_taiwan.json', 'w')); print(len(g), g[0], g[-1])
    elif mode == 'test':
        g = json.load(open('granules_l4a_taiwan.json')); D = defs(); json.dump(D, open('definitions.json', 'w'))
        for x in g[100:103]:
            res, sc, n, dt_, body = fetch(x['name']); a, m = decode(body, D)
            print(res, sc, n, dt_, 's shots', len(a), [t for c, t in m][-2:]);
            if len(a): print(a[:2], 'agbd range', float(a['agbd'].min()), float(a['agbd'].max()), 'flags', np.unique(a['flags']))
    elif mode == 'run':
        out = sys.argv[2]; os.makedirs(out + '/responses', exist_ok=True)
        g = json.load(open('granules_l4a_taiwan.json')); log = open(out + '/fetch_log.jsonl', 'a')
        done = {f[:-4] for f in os.listdir(out + '/responses') if f.endswith('.bin')}
        todo = [x['name'] for x in g if (x['name'].split('.', 1)[1] if x['name'].startswith('GEDI_L4A') else x['name']) not in done]
        print('todo', len(todo), 'done', len(done), flush=True)
        DEADLINE = time.time() + float(sys.argv[4]) if len(sys.argv) > 4 else 1e18
        def job(nm):
            if time.time() > DEADLINE: return None
            for attempt in range(2):
                try:
                    res, sc, n, dt_, body = fetch(nm)
                    if sc == 200:
                        tmp = out + '/responses/' + res + '.part'; open(tmp, 'wb').write(body); os.replace(tmp, out + '/responses/' + res + '.bin')
                        return {'res': res, 'status': sc, 'bytes': n, 'sec': dt_, 'sha256': hashlib.sha256(body).hexdigest()}
                    err = {'res': res, 'status': sc, 'bytes': n, 'sec': dt_}
                except Exception as e:
                    err = {'res': nm, 'error': str(e)[:200]}
                time.sleep(3)
            return err
        with ThreadPoolExecutor(int(sys.argv[3]) if len(sys.argv) > 3 else 4) as ex:
            for k, f in enumerate(as_completed([ex.submit(job, n) for n in todo])):
                if f.result() is None: continue
                log.write(json.dumps(f.result()) + '\n'); log.flush()
                if k % 25 == 0: print(k, f.result(), flush=True)
        print('RUN DONE', flush=True)
    elif mode == 'decode':
        out = sys.argv[2]; D = json.load(open('definitions.json')); parts = []; bad = 0; nres = 0
        for f in sorted(os.listdir(out + '/responses')):
            if not f.endswith('.bin'): continue
            a, m = decode(open(out + '/responses/' + f, 'rb').read(), D); nres += 1
            if not any('Successfully completed' in t for c, t in m): bad += 1
            if len(a): parts.append(a)
        A = np.concatenate(parts); print('granules', nres, 'without completion message', bad, 'shots', len(A))
        np.savez_compressed(sys.argv[3], **{k: A[k] for k in A.dtype.names})
