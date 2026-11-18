# -*- coding: utf-8 -*-
"""s09: supporting tables and statistics (recovery bins, scar-size gradient, official annual landslide maps, GEDI checks, catalogue summary)."""
import json, numpy as np, pandas as pd
D = 'v2'; TG = {1: 'BL', 2: 'CF', 4: 'MX', 3: 'BB', 5: 'BB', 6: 'BB'}; GS = ['BL', 'CF', 'MX', 'BB']; OUT = {}
P = json.load(open(D + '/parameters.json')); S = json.load(open(D + '/ledger_summary.json'))
A = pd.read_parquet(D + '/cells.parquet')
A['tg'] = A.type_nb.map(TG); A['eb'] = np.digitize(A.elev.fillna(0), [500, 1500, 2500]); A['tcb'] = np.digitize(A.tc2000.fillna(0), [25.001, 50.001, 75.001, 90.001])
A['blk'] = (A.x // 10000).astype(int) * 1000 + (A.y // 10000).astype(int); A['iy'] = np.floor(A.t_chm); A['sz6'] = np.digitize(A.poly_ha.fillna(0), [1, 3, 10, 30, 100]); A['sz'] = np.digitize(A.poly_ha.fillna(0), [3, 30])
R = A[A.is_ref & A.chm.notna() & A.t_chm.notna()]
r1 = R.groupby(['blk', 'iy', 'tg']).chm.agg(href_t='mean', n_t='size').reset_index(); r2 = R.groupby(['blk', 'iy']).chm.agg(href_a='mean', n_a='size').reset_index()
A = A.merge(r1, on=['blk', 'iy', 'tg'], how='left').merge(r2, on=['blk', 'iy'], how='left'); A['href'] = np.where(A.n_t >= 30, A.href_t, np.where(A.n_a >= 30, A.href_a, np.nan)); A['xr'] = A.chm / A.href
L = A[~A.is_ref]; core = L[(L.a >= 0.8) & L.xr.notna() & L.tg.notna()]
post = core[(core.t_chm > core.t_last_after) & (core.nyears == 1)].copy(); post['age'] = post.t_chm - post.t_last; pre = core[core.t_chm < core.t_first_before].copy(); F = post[post.tcb >= 3]
# --- recovery bins (forest-origin, single failure)
bins = [0, 1, 2, 4, 6, 8, 10, 12, 16]; F = F.assign(ab=pd.cut(F.age, bins))
Fb = F[~((F.tg == 'BB') & (F.sz == 2))]      # bamboo: scars < 30 ha, as in the rate regression
OUT['bb_large'] = {'cells': int(((F.tg == 'BB') & (F.sz == 2)).sum()), 'distinct_polygon_sizes': int(F[(F.tg == 'BB') & (F.sz == 2)].poly_ha.round(3).nunique()), 'x_age_lt5': float(F[(F.tg == 'BB') & (F.sz == 2) & (F.age < 5)].xr.mean()), 'n_age_lt5': int(((F.tg == 'BB') & (F.sz == 2) & (F.age < 5)).sum())}
rb = Fb.groupby(['tg', 'ab'], observed=True).agg(age=('age', 'mean'), x=('xr', 'mean'), h=('chm', 'mean'), href=('href', 'mean'), lt2=('chm', lambda v: (v < 2).mean()), n=('xr', 'size')).reset_index(); rb['ab'] = rb.ab.astype(str); rb.to_csv(D + '/recovery_bins.csv', index=False)
T3 = F[F.tg != 'BB']
rs = T3.assign(szl=T3.sz6).groupby(['szl', 'ab'], observed=True).agg(age=('age', 'mean'), x=('xr', 'mean'), n=('xr', 'size')).reset_index(); rs['ab'] = rs.ab.astype(str); rs.to_csv(D + '/recovery_bins_by_size.csv', index=False)
st = {}
for g in GS:
    d = F[F.tg == g]; y = d[d.age < 2]; o = d[(d.age >= 8) & (d.age < 12)]; q = pre[(pre.tg == g) & (pre.tcb >= 3)]
    st[g] = {'h_lt2yr': float(y.chm.mean()), 'h_8_12': float(o.chm.mean()), 'x_8_12': float(o.xr.mean()), 'lt2m_8_12': float((o.chm < 2).mean()), 'n_8_12': int(len(o)), 'age_8_12': float(o.age.mean()), 'h_pre_tc75': float(q.chm.mean()), 'x_pre_tc75': float(q.xr.mean()), 'href_pre': float(q.href.mean()), 'n_pre_tc75': int(len(q)),
             'rec_frac_8_12': float((o.xr.mean() - P['x0'][g]) / (q.xr.mean() - P['x0'][g])), 'years_to_pre_tc75_90': float((1 - P['x0'][g]) * P['rho'][g][3] / P['k'][g]), 'years_to_pre_tc90': float((1 - P['x0'][g]) * P['rho'][g][4] / P['k'][g]), 'years_to_reference': float((1 - P['x0'][g]) / P['k'][g])}
OUT['recovery_stats'] = st; OUT['lt2m_8_12_all'] = float((F[(F.age >= 8) & (F.age < 12)].chm < 2).mean()); OUT['n_forest_origin_post'] = int(len(F)); OUT['n_post_single'] = int(len(post)); OUT['n_pre'] = int(len(pre))
OUT['ref_height_by_image_year'] = R.groupby('iy').chm.mean().round(2).to_dict(); OUT['chm_year_share'] = (core.iy.value_counts(normalize=True).sort_index().round(4)).to_dict(); OUT['t_chm_median'] = float(core.t_chm.median())
# pre-failure relative biomass by year of first failure (closed-canopy cells, all tree types pooled and broadleaf)
x0m = pre.tg.map(P['x0']); pre = pre.assign(rho=((pre.xr - x0m) / (1 - x0m)).clip(lower=0))
py = pre[pre.tcb >= 3].groupby('first').agg(n=('rho', 'size'), rho=('rho', 'mean'), h=('chm', 'mean'), href=('href', 'mean')).reset_index(); py = py[py.n >= 200]; py.to_csv(D + '/rho_by_failure_year.csv', index=False)
# --- scar-size gradient: regression rate, relative height at 8-12 yr, official bare share
def fe(d, keys):
    g = d.groupby(keys); a = d.age - g.age.transform('mean'); h = d.xr - g.xr.transform('mean'); return float((a * h).sum() / (a * a).sum())
I1 = L[(L.a >= 0.8) & (L.a3 >= 0.6) & (L.nyears == 1) & L.tg.notna()].copy(); IF = I1[I1.tcb >= 3]; YRS = list(range(2017, 2025))
rng = np.random.default_rng(7); blocks = T3.blk.unique(); gb = {b: d for b, d in T3.groupby('blk')}
boot = []
for it in range(200):
    d = pd.concat([gb[b] for b in rng.choice(blocks, len(blocks))]); boot.append([fe(d[d.sz6 == q], ['eb', 'tg']) if (d.sz6 == q).sum() > 300 else np.nan for q in range(6)])
boot = np.array(boot); rows = []
for q, lab in enumerate(['<1', '1-3', '3-10', '10-30', '30-100', '>=100']):
    d = T3[T3.sz6 == q]; e = IF[(IF.sz6 == q) & (IF.tg != 'BB') & (IF['last'] <= 2016)]; e9 = e[e['last'] == 2009]
    rows.append({'size_ha': lab, 'n_cells': len(d), 'k': fe(d, ['eb', 'tg']), 'k_lo': np.nanpercentile(boot[:, q], 2.5), 'k_hi': np.nanpercentile(boot[:, q], 97.5), 'x_8_12': d[(d.age >= 8) & (d.age < 12)].xr.mean(), 'n_8_12': int(((d.age >= 8) & (d.age < 12)).sum()),
                 'bare2017': (e.s2017 >= 0.5).mean(), 'bare2024': (e.s2024 >= 0.5).mean(), 'n_map': len(e), 'bare2017_2009': (e9.s2017 >= 0.5).mean(), 'bare2024_2009': (e9.s2024 >= 0.5).mean(), 'n_map_2009': len(e9)})
pd.DataFrame(rows).to_csv(D + '/size_gradient.csv', index=False)
# --- official annual landslide maps: persistence of bare state on forest-origin scars
per = IF[IF['last'] <= 2016].groupby('last').apply(lambda d: pd.Series({'n': len(d), **{str(y): (d[f's{y}'] >= 0.5).mean() for y in YRS}})).reset_index(); per.to_csv(D + '/persistence_by_failure_year.csv', index=False)
e9 = IF[IF['last'] == 2009]
pt = e9.groupby('tg').apply(lambda d: pd.Series({'n': len(d), **{str(y): (d[f's{y}'] >= 0.5).mean() for y in YRS}})).reset_index(); pt.to_csv(D + '/persistence_2009_by_type.csv', index=False)
M9 = np.column_stack([(e9[f's{y}'] >= 0.5).values for y in YRS])
OUT['persistence_2009'] = {'n': int(len(e9)), 'share': {str(y): float(M9[:, i].mean()) for i, y in enumerate(YRS)}, 'ever': float(M9.any(1).mean()), 'all8': float(M9.all(1).mean()), 'never': float((~M9.any(1)).mean()),
                           'bare_to_not': [float((M9[:, i] & ~M9[:, i + 1]).sum() / M9[:, i].sum()) for i in range(7)], 'not_to_bare': [float((~M9[:, i] & M9[:, i + 1]).sum() / (~M9[:, i]).sum()) for i in range(7)],
                           'slope_pp_per_yr': float(np.polyfit(YRS, 100 * M9.mean(0), 1)[0])}
OUT['ref_mapped_share'] = {str(y): float((A[A.is_ref][f's{y}'] >= 0.5).mean()) for y in YRS}
# CHM by official state in the image year
rows = []
cs = L[(L.a >= 0.8) & L.chm.notna() & L.tg.notna() & (L.tcb >= 3)]
for y in (2017, 2018, 2019):
    d = cs[(cs.t_chm >= y) & (cs.t_chm < y + 1) & (cs.t_last_after < y)]
    for lab, sel in [('mapped as landslide', d[f's{y}'] >= 0.5), ('not mapped', d[f's{y}'] < 0.1)]:
        rows.append({'year': y, 'state': lab, 'n': int(sel.sum()), 'chm_mean': d[sel].chm.mean(), 'chm_median': d[sel].chm.median(), 'x_mean': d[sel].xr.mean(), 'share_lt2m': (d[sel].chm < 2).mean(), 'age': (d[sel].t_chm - d[sel].t_last).mean()})
pd.DataFrame(rows).to_csv(D + '/chm_by_official_state.csv', index=False)
d = cs[(cs.t_chm >= 2017) & (cs.t_chm < 2020) & (cs.t_last_after < 2017)]; d = d.assign(st=np.select([sum((d.iy == y) * (d[f's{y}'] >= 0.5) for y in (2017, 2018, 2019)) > 0, sum((d.iy == y) * (d[f's{y}'] < 0.1) for y in (2017, 2018, 2019)) > 0], ['bare', 'not'], 'mixed'))
OUT['chm_by_state_pooled'] = {k: {'n': int(len(v)), 'chm': float(v.chm.mean()), 'lt2': float((v.chm < 2).mean()), 'x': float(v.xr.mean())} for k, v in d.groupby('st')}
# agreement of the event catalogue with the annual maps (single-failure interior cells first catalogued 2017-2024)
rows = []
for y0 in range(2017, 2025):
    d = I1[I1['first'] == y0]; r = {'first': y0, 'n': len(d)}
    for kk in (-1, 0, 1):
        y = y0 + kk; r[f'map{kk:+d}'] = (d[f's{y}'] >= 0.5).mean() if 2017 <= y <= 2024 else np.nan
    rows.append(r)
ag = pd.DataFrame(rows); ag['same_or_next'] = ag[['map+0', 'map+1']].max(axis=1); ag.to_csv(D + '/catalogue_vs_annual_map.csv', index=False)
OUT['agreement_weighted'] = float(np.average(ag.same_or_next, weights=ag.n)); OUT['mapped_before_weighted'] = float(np.average(ag['map-1'].fillna(0)[1:], weights=ag.n[1:]))
# GFC loss-year timing on single-event cores
c1 = L[(L.a >= 0.8) & (L.nyears == 1) & (L.tc2000 >= 50)]; dd = (2000 + c1.lossyear.where(c1.lossyear > 0)) - c1['first']
OUT['gfc'] = {'n': int(len(c1)), 'flagged': float(c1.lossyear.gt(0).mean()), 'same_year_given_flag': float((dd == 0).sum() / dd.notna().sum()), 'within1_given_flag': float(dd.between(-1, 1).sum() / dd.notna().sum()), 'ref_flagged': float((A[A.is_ref & (A.tc2000 >= 50)].lossyear > 0).mean())}
# --- GEDI L4A checks
G = pd.read_parquet(D + '/gedi_l4a_joined.parquet'); G['tg'] = G.type.map(TG); G['sb'] = np.digitize(G.slope.fillna(0), [10, 20, 30, 40])
Rg = G[G.tg.notna() & ~G.ls5 & G.chm.notna()]; bl = Rg[Rg.tg == 'BL'].assign(hb=pd.cut(Rg[Rg.tg == 'BL'].chm, [0, 2, 4, 6, 8, 10, 12, 15, 18, 22, 30]))
gt = bl.groupby(['hb', 'sb'], observed=True).agg(agbd=('agbd', 'mean'), h=('chm', 'mean'), n=('agbd', 'size')).reset_index(); gt['hb'] = gt.hb.astype(str); gt.to_csv(D + '/gedi_agbd_by_height_slope_BL.csv', index=False)
g0 = bl[bl.sb == 0]; lin = np.polyfit(g0.chm, g0.agbd, 1); bm = g0.groupby('hb', observed=True).agg(h=('chm', 'mean'), a=('agbd', 'mean')); lb = np.polyfit(bm.h, bm.a, 1)
Sg = G[(G.a >= 0.8) & (G.a3 >= 0.5) & G.chm.notna()]; pg = Sg[(Sg.n_before > 0) & (Sg.n_amb == 0) & (Sg.n_after.fillna(0) == 0)].copy(); pg['age'] = pg.t - pg.t_last; pg['s3'] = np.digitize(pg.slope.fillna(0), [20, 30])
gs = pg.assign(ab=pd.cut(pg.age, [0, 2, 5, 10, 13, 16, 22])).groupby(['ab', 's3'], observed=True).agg(agbd=('agbd', 'mean'), chm=('chm', 'mean'), n=('agbd', 'size')).reset_index(); gs['ab'] = gs.ab.astype(str); gs.to_csv(D + '/gedi_scar_agbd_by_age_slope.csv', index=False)
prg = Sg[(Sg.n_before == 0) & (Sg.n_amb == 0) & (Sg.n_after > 0)]
OUT['gedi'] = {'shots': int(len(G)), 'ref_shots': int(len(Rg)), 'ref_BL': int(len(bl)), 'scar_post': int(len(pg)), 'scar_pre': int(len(prg)), 'lin_slope_lt10deg': [float(lin[0]), float(lin[1])], 'lin_binned_lt10deg': [float(lb[0]), float(lb[1])], 'r_lt10': float(np.corrcoef(g0.chm, g0.agbd)[0, 1]), 'n_lt10': int(len(g0)),
               'agbd_h0_2_by_slope': bl[bl.chm <= 2].groupby('sb').agbd.mean().round(1).tolist(), 'agbd_h12_15_by_slope': bl[(bl.chm > 12) & (bl.chm <= 15)].groupby('sb').agbd.mean().round(1).tolist(), 'scar_slope_median': float(pg.slope.median()), 'scar_lt20_share': float((pg.slope < 20).mean()),
               'scar_agbd_lt2yr': float(pg[pg.age < 2].agbd.mean()), 'scar_agbd_10_16': float(pg[(pg.age >= 10) & (pg.age < 16)].agbd.mean()), 'scar_chm_10_16': float(pg[(pg.age >= 10) & (pg.age < 16)].chm.mean()), 'pre_agbd': float(prg.agbd.mean()), 'pre_chm': float(prg.chm.mean()), 't_min': float(G.t.min()), 't_max': float(G.t.max())}
# --- catalogue summary and unit conversions
from config import WORK
rec = pd.read_csv(WORK + '/q01/q01_records.csv')
cat = rec.groupby('catalog_year').agg(polygons=('fid', 'size'), area_ha=('area_best_ha', 'sum'), raster_ha=('area_raster_ha', 'sum'), median_ha=('area_best_ha', 'median'), max_ha=('area_best_ha', 'max')).reset_index(); cat.to_csv(D + '/catalogue_by_year.csv', index=False)
OUT['catalogue'] = {'polygons': int(len(rec)), 'area_ha': float(rec.area_best_ha.sum()), 'raster_ha': float(rec.area_raster_ha.sum()), 'median_ha': float(rec.area_best_ha.median()), 'share_ge10ha_area': float(rec[rec.area_best_ha >= 10].area_best_ha.sum() / rec.area_best_ha.sum()), 'n_ge10ha': int((rec.area_best_ha >= 10).sum())}
Bk = pd.read_csv(D + '/ledger_breakdown.csv'); t9 = Bk[(Bk.dim == 'type') & (Bk.year == 2009)].set_index('class').loss_Mt; FAC = {'BL': 1.24 * 0.4691 * 44 / 12, 'CF': 1.22 * 0.4821 * 44 / 12, 'MX': 1.23 * 0.4756 * 44 / 12, 'BB': 1.46 * 0.4732 * 44 / 12}
OUT['agb_dm_2009_Tg'] = float(sum(t9.get(g, 0) / FAC[g] for g in GS)); OUT['agb_dm_2009_variants_Tg'] = {k: float(v['loss_2009'] / S['loss_2009'] * OUT['agb_dm_2009_Tg']) for k, v in S['variants'].items()}
OUT['carbon_2009_TgC'] = float(S['loss_2009'] * 12 / 44); OUT['volume_equiv_2009_Mm3'] = float(sum(t9.get(g, 0) / r for g, r in [('BL', 1.9622), ('CF', 1.0999), ('MX', 1.5310)]))
st_ = pd.read_csv('v1/nfi4_strata_type_elev.csv'); st_.to_csv(D + '/nfi4_strata_type_elev.csv', index=False)
v1 = json.load(open('v1/ledger_summary.json')); OUT['hv'] = v1.get('hv'); OUT['v1'] = {'loss_total': v1['loss_total'], 'gain_total': v1['gain_total']}
json.dump(OUT, open(D + '/support_stats.json', 'w'), indent=1)
pd.set_option('display.width', 250)
print(json.dumps(OUT['recovery_stats'], indent=0)[:2200]); print(pd.read_csv(D + '/size_gradient.csv').round(4).to_string(index=False)); print(per.round(3).to_string(index=False)); print(pt.round(3).to_string(index=False))
print({k: v for k, v in OUT.items() if k in ('persistence_2009', 'chm_by_state_pooled', 'agreement_weighted', 'mapped_before_weighted', 'gfc', 'gedi', 'catalogue', 'agb_dm_2009_Tg', 'carbon_2009_TgC', 'volume_equiv_2009_Mm3', 'lt2m_8_12_all', 't_chm_median')})
print(ag.round(3).to_string(index=False)); print(pd.read_csv(D + '/chm_by_official_state.csv').round(2).to_string(index=False))
