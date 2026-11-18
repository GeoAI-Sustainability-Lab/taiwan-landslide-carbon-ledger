# -*- coding: utf-8 -*-
"""s08: national living-biomass ledger of catalogued landslides, Taiwan 2004-2025.
Pre-failure stock = official stand stock x relative pre-failure canopy stature of failing cells; recovery from locally
normalised canopy height on scars (CHMv2); block bootstrap + official conversion-factor uncertainty; structural variants."""
import os, sys, json, numpy as np, pandas as pd
D = 'v2'; NB = int(os.environ.get('NBOOT', 300)); rng = np.random.default_rng(20261111)
TG = {1: 'BL', 2: 'CF', 4: 'MX', 3: 'BB', 5: 'BB', 6: 'BB'}; GS = ['BL', 'CF', 'MX', 'BB']
RATIO = {1: 1.9622, 2: 1.0999, 4: 1.5310, 5: 2.5256, 6: 2.1150}; BAMBOO_CO2 = 176.5
U95 = {'BL': 0.1224, 'CF': 0.0686, 'MX': 0.0844, 'BB': 0.1224}      # 95% half-width of BCEF x (1+R) x CF, from NIR 2026 Table 6.2.9
A = pd.read_parquet(D + '/cells.parquet'); CY = pd.read_parquet('v1/cell_year_events.parquet')
ST = pd.read_csv('v1/nfi4_strata_type_elev.csv').set_index(['type', 'band']).vol_mean.to_dict()
A['tg'] = A.type_nb.map(TG); A['eb'] = np.digitize(A.elev.fillna(0), [500, 1500, 2500]); A['tcb'] = np.digitize(A.tc2000.fillna(0), [25.001, 50.001, 75.001, 90.001])
A['blk'] = (A.x // 10000).astype(int) * 1000 + (A.y // 10000).astype(int); A['iy'] = np.floor(A.t_chm)
A['sz'] = np.digitize(A.poly_ha.fillna(0), [3, 30])                 # scar size class: <3, 3-30, >=30 ha
A['sz6'] = np.digitize(A.poly_ha.fillna(0), [1, 3, 10, 30, 100])
# ---- local reference height: never-catalogued forest of the same type group, 10 km block and image year
R = A[A.is_ref & A.chm.notna() & A.t_chm.notna()]
r1 = R.groupby(['blk', 'iy', 'tg']).chm.agg(href_t='mean', n_t='size').reset_index(); r2 = R.groupby(['blk', 'iy']).chm.agg(href_a='mean', n_a='size').reset_index()
A = A.merge(r1, on=['blk', 'iy', 'tg'], how='left').merge(r2, on=['blk', 'iy'], how='left')
A['href'] = np.where(A.n_t >= 30, A.href_t, np.where(A.n_a >= 30, A.href_a, np.nan)); A['xr'] = A.chm / A.href
L = A[~A.is_ref].reset_index(drop=True); Rf = A[A.is_ref & A.xr.notna()].copy()
def samples(amin=0.8, interior=None, single=True):
    core = L[(L.a >= amin) & L.xr.notna() & L.tg.notna()]
    if interior is not None: core = core[core.a3 >= interior]
    post = core[core.t_chm > core.t_last_after].copy(); post['age'] = post.t_chm - post.t_last
    if single: post = post[post.nyears == 1].copy()      # baseline and recovery from cells with one catalogued failure year
    pre = core[core.t_chm < core.t_first_before].copy()
    return post, pre
post, pre = samples()
def fe(d, keys):
    g = d.groupby(keys); a = d.age - g.age.transform('mean'); h = d.xr - g.xr.transform('mean'); return float((a * h).sum() / (a * a).sum())
def params(po, pr, level=False):
    """x0: relative height on scars imaged < 1 yr after failure; k: annual gain of relative height on forest-origin scars (regression on years since
    failure with elevation-band and scar-size-class effects); rho[g, c]: relative pre-failure biomass of failing cells by type group and 2000 tree-cover class"""
    x0 = np.zeros(4); k = np.zeros(4); rho = np.zeros((4, 5)); x0p = po[po.age < 1].xr.mean()
    F = po[po.tcb >= 3]
    for j, g in enumerate(GS):
        y = po[(po.tg == g) & (po.age < 1)].xr; x0[j] = y.mean() if len(y) >= 30 else x0p
        d = F[F.tg == g]
        if g == 'BB': d = d[d.sz < 2]
        if level: f = d[(d.age >= 8) & (d.age < 12)]; k[j] = (f.xr.mean() - x0[j]) / f.age.mean() if len(f) >= 30 else np.nan
        else: k[j] = fe(d, ['eb', 'sz6']) if len(d) >= 500 else np.nan
    k = np.where(np.isfinite(k), k, np.nanmean(k))
    for c in range(5):
        pc = pr[pr.tcb == c]; pool = np.mean([((pc[pc.tg == g].xr - x0[j]) / (1 - x0[j])).clip(lower=0).mean() for j, g in enumerate(GS) if (pc.tg == g).sum() >= 50] or [0.0])
        for j, g in enumerate(GS):
            y = pc[pc.tg == g].xr; rho[j, c] = ((y - x0[j]) / (1 - x0[j])).clip(lower=0).mean() if len(y) >= 50 else pool
    return x0, np.clip(k, 0, None), rho
x0, k, rho = params(post, pre)
rho_ref = np.array([[((Rf[(Rf.tg == g) & (Rf.tcb == c)].xr - x0[j]) / (1 - x0[j])).clip(lower=0).mean() for c in range(5)] for j, g in enumerate(GS)])
# ---- cell-level base stock (official): stand polygon stock where mapped, else type x elevation stratum mean
tc = L.type_cell.values; tnb = L.type_nb.values; vol = L.vol.values; eb = L.eb.values
def strat(types):
    s = np.zeros(len(L))
    for t in (1, 2, 4, 5, 6):
        v = np.array([ST.get((t, b), np.nan) for b in range(4)]); v = np.where(np.isfinite(v), v, np.nanmean(v)); s[types == t] = (v[eb] * RATIO[t])[types == t]
    s[types == 3] = BAMBOO_CO2; return s
typed = (tc >= 1) & (tc <= 6)
S_poly = np.zeros(len(L))
for t in (1, 2, 4, 5, 6): m = (tc == t) & (vol > 0); S_poly[m] = vol[m] * RATIO[t]
S_base = np.where(typed & (S_poly > 0), S_poly, np.where(typed, strat(tc), strat(tnb))); S_strat = np.where(typed, strat(tc), strat(tnb))
hasg = L.tg.notna().values; gi = np.array([GS.index(g) if isinstance(g, str) else 0 for g in L.tg.values]); ci = L.tcb.values; first = L['first'].values
area = L.a.values * 0.09
pos = pd.Series(np.arange(len(L)), index=L.cell.values); CY['i'] = pos.reindex(CY.cell.values).values; CY['p'] = (CY.frac.values / L.a.values[CY.i.values]).clip(0, 1)
years = np.arange(2004, 2026); evy = {y: (d.i.values, d.p.values) for y, d in CY.groupby('year')}
def run(S, rate, detail=False):
    f = np.ones(len(L)); touched = np.zeros(len(L), bool); loss = np.zeros(22); gain = np.zeros(22); out = {}
    if detail: lc = np.zeros((22, len(L)), np.float32) if False else None; lossc = np.zeros(len(L)); gainc = np.zeros(len(L)); lossy = {}
    for kk, y in enumerate(years):
        if touched.any():
            fn = np.minimum(1.0, f + rate); g_ = np.where(touched, (fn - f) * S * area, 0.0); gain[kk] = g_.sum(); f = np.where(touched, fn, f)
            if detail: gainc += g_
        if y in evy:
            i, p = evy[y]; l_ = p * f[i] * S[i] * area[i]; loss[kk] = l_.sum()
            if detail: lossc[i] += l_; lossy[y] = (i, l_)
            f[i] = f[i] * (1 - p); touched[i] = True
    deficit = ((1 - f) * S * area)[touched].sum()
    if detail: return loss / 1e6, gain / 1e6, deficit / 1e6, lossc, gainc, lossy, f
    return loss / 1e6, gain / 1e6, deficit / 1e6
def build(x0_, k_, rho_, rho_early=None, conv=None, Sb=S_base, kcell=None, x0cell=None):
    rc = np.where(hasg, rho_[gi, ci], 0.0)
    if rho_early is not None: rc = np.where(hasg & (first <= 2010), rho_early[gi, ci], rc)
    S = Sb * rc * (1.0 if conv is None else conv[gi])
    kc = k_[gi] if kcell is None else kcell; xc = x0_[gi] if x0cell is None else x0cell
    rate = np.where(rc > 0, kc / ((1 - xc) * np.maximum(rc, 1e-6)), 0.0)
    return S, rate
S_c, rate_c = build(x0, k, rho)
loss_c, gain_c, def_c, lossc, gainc, lossy, f_end = run(S_c, rate_c, detail=True)
# ---- bootstrap: 10 km blocks (stature, baseline, recovery) x conversion-factor uncertainty
pb = {b: d for b, d in post.groupby('blk')}; qb = {b: d for b, d in pre.groupby('blk')}; blocks = np.array(sorted(set(pb) | set(qb)))
bl, bg, bd, bx0, bk, brho = [], [], [], [], [], []
for it in range(NB):
    bs = rng.choice(blocks, len(blocks)); po = pd.concat([pb[b] for b in bs if b in pb]); pr = pd.concat([qb[b] for b in bs if b in qb])
    x0b, kb, rb = params(po, pr); conv = np.array([1 + rng.normal(0, U95[g] / 1.96) for g in GS]); conv[3] = conv[0]
    S, rt = build(x0b, kb, rb, conv=conv); l, g, d = run(S, rt); bl.append(l); bg.append(g); bd.append(d); bx0.append(x0b); bk.append(kb); brho.append(rb)
bl = np.array(bl); bg = np.array(bg); bd = np.array(bd); bx0 = np.array(bx0); bk = np.array(bk); brho = np.array(brho)
ci95 = lambda v, ax=0: np.percentile(v, [2.5, 97.5], axis=ax)
# ---- structural variants (central parameters unless stated)
F = post[post.tcb >= 3]
kspec = {}
for nm, keys in [('elev', ['eb']), ('elev+size', ['eb', 'sz6']), ('elev+image year', ['eb', 'iy']), ('elev+image year+size', ['eb', 'iy', 'sz6']), ('elev+failure year', ['eb', 'last']), ('elev+failure year+size', ['eb', 'last', 'sz6']), ('elev+block+size', ['eb', 'blk', 'sz6'])]:
    kspec[nm] = {g: fe(F[F.tg == g], keys) for g in GS[:3]}; kspec[nm]['ALL'] = fe(F, keys + ['tg'])
V = {}
def rec(name, S, rate):
    l, g, d = run(S, rate); V[name] = {'loss_total': float(l.sum()), 'loss_2009': float(l[5]), 'gain_total': float(g.sum()), 'gain_2025': float(g[-1]), 'deficit_2025': float(d), 'offset_pct': float(100 * g.sum() / l.sum())}
rec('central', S_c, rate_c)
Su = S_base * np.where(hasg & (L.tc2000.fillna(0).values >= 30), 1.0, 0.0); rec('unadjusted stand stock (tree cover 2000 >= 30%)', Su, np.where(Su > 0, k[gi] / (1 - x0[gi]), 0.0))
rec('reference stature for failures up to 2010', *build(x0, k, rho, rho_early=rho_ref))
rec('stratum-mean stock for all cells', *build(x0, k, rho, Sb=S_strat))
klo = np.array([kspec['elev+failure year'][g] for g in GS[:3]] + [k[3]]); khi = np.array([kspec['elev+image year+size'][g] for g in GS[:3]] + [k[3]])
rec('recovery rate: within failure-year regression', *build(x0, np.clip(klo, 0, None), rho)); rec('recovery rate: within image-year and size regression', *build(x0, khi, rho))
a_, b_, c_ = params(post, pre, level=True); rec('recovery rate: level at 8-12 yr minus baseline', *build(a_, b_, c_)); V['recovery rate: level at 8-12 yr minus baseline']['k_BL'] = float(b_[0])
# size-stratified baseline and recovery (three scar-size classes), pooled over tree types
x0s = np.zeros(3); ks = np.zeros(3)
for s in range(3):
    d = post[(post.sz == s) & (post.tg != 'BB')]; x0s[s] = d[d.age < 2].xr.mean(); ks[s] = max(fe(d[d.tcb >= 3], ['eb', 'tg']), 0)
szc = L.sz.values; rec('scar-size-specific recovery', *build(x0, k, rho, kcell=np.where(gi == 3, k[3], ks[szc]), x0cell=np.where(gi == 3, x0[3], x0s[szc])))
for nm, kw in [('parameters from cells with landslide fraction >= 0.5', dict(amin=0.5)), ('parameters from cells with landslide fraction >= 0.95', dict(amin=0.95)), ('parameters from scar interiors (3x3 mean fraction >= 0.8)', dict(interior=0.8))]:
    po, pr = samples(**kw); a_, b_, c_ = params(po, pr); rec(nm, *build(a_, b_, c_)); V[nm]['k_BL'] = float(b_[0]); V[nm]['rho_BL_tc90'] = float(c_[0, 4])
po_, pr_ = samples(single=False); a_, b_, c_ = params(po_, pr_); rec('recovery parameters including cells with repeated failures', *build(a_, b_, c_)); V['recovery parameters including cells with repeated failures']['k_BL'] = float(b_[0])
# quadratic biomass-height scaling (b = 2)
def params_b2(po, pr):
    r2 = np.zeros((4, 5)); nrm = np.zeros(4)
    for j, g in enumerate(GS):
        nrm[j] = (((Rf[Rf.tg == g].xr - x0[j]) / (1 - x0[j])).clip(lower=0) ** 2).mean()
        for c in range(5):
            y = pr[(pr.tg == g) & (pr.tcb == c)].xr; r2[j, c] = ((((y - x0[j]) / (1 - x0[j])).clip(lower=0) ** 2).mean() / nrm[j]) if len(y) >= 50 else np.nan
    for c in range(5): r2[:, c] = np.where(np.isfinite(r2[:, c]), r2[:, c], np.nanmean(r2[:, c]))
    return r2, nrm
rho2, nrm2 = params_b2(post, pre)
def run_b2():
    rc = np.where(hasg, rho2[gi, ci], 0.0); S = S_base * rc; xe = np.sqrt(np.maximum(rc * nrm2[gi], 1e-9))      # effective pre-failure relative height above baseline
    rate = np.where(rc > 0, k[gi] / ((1 - x0[gi]) * xe), 0.0)
    ph = np.ones(len(L)); touched = np.zeros(len(L), bool); loss = np.zeros(22); gain = np.zeros(22)
    for kk, y in enumerate(years):
        if touched.any():
            pn = np.minimum(1.0, ph + rate); g_ = np.where(touched, (pn ** 2 - ph ** 2) * S * area, 0.0); gain[kk] = g_.sum(); ph = np.where(touched, pn, ph)
        if y in evy:
            i, p = evy[y]; loss[kk] = (p * ph[i] ** 2 * S[i] * area[i]).sum(); ph[i] = np.sqrt(ph[i] ** 2 * (1 - p)); touched[i] = True
    d = ((1 - ph ** 2) * S * area)[touched].sum(); return loss / 1e6, gain / 1e6, d / 1e6
l, g, d = run_b2(); V['quadratic biomass-height scaling'] = {'loss_total': float(l.sum()), 'loss_2009': float(l[5]), 'gain_total': float(g.sum()), 'gain_2025': float(g[-1]), 'deficit_2025': float(d), 'offset_pct': float(100 * g.sum() / l.sum())}
# forward run without new failures (constant recovery rates)
fw = f_end.copy(); tch = (lossc > 0) | (f_end < 1); proj = {}
for y in range(2026, 2126):
    fw = np.where(tch, np.minimum(1.0, fw + rate_c), fw)
    if y in (2030, 2040, 2050, 2075, 2100, 2125): proj[y] = float(((1 - fw) * S_c * area)[tch].sum() / 1e6)
    if 'half' not in proj and ((1 - fw) * S_c * area)[tch].sum() / 1e6 <= def_c / 2: proj['half'] = y
# ---- official inventory series (NIR 2026, chapter 6): forest-land net removals and disturbance losses, kt CO2
NIR_REM = dict(zip(range(2004, 2025), [22542, 22290, 22259, 22074, 22082, 19388, 21889, 21947, 21960, 21974, 21886, 21900, 21926, 21961, 21980, 21915, 21899, 21843, 21826, 21717, 21757]))
NIR_DIST = dict(zip(range(2004, 2025), [61, 199, 8, 71, 8, 2535, 8, 0.5, 0.2, 0.7, 0.8, 1, 5, 0.7, 4, 1, 1, 27, 29, 81, 0.5]))
NIR_AREA = dict(zip(range(2010, 2025), [12109, 2821, 3197, 2104, 2408, 1080, 834, 1500, 804, 1521, 890, 120, 524, 620, 757]))
# ---- tables
T = pd.DataFrame({'year': years, 'loss': loss_c, 'loss_lo': ci95(bl)[0], 'loss_hi': ci95(bl)[1], 'gain': gain_c, 'gain_lo': ci95(bg)[0], 'gain_hi': ci95(bg)[1]})
T['net'] = T.gain - T.loss; T['nir_removal'] = [NIR_REM.get(y, np.nan) / 1e3 for y in years]; T['nir_disturbance'] = [NIR_DIST.get(y, np.nan) / 1e3 for y in years]; T['loss_pct_removal'] = 100 * T.loss / T.nir_removal
ar = pd.DataFrame({'year': CY.year.values, 'ha': CY.frac.values * 0.09, 'i': CY.i.values}); ar['tcb'] = ci[ar.i.values]; ar['hasg'] = hasg[ar.i.values]
T['area_ha'] = ar.groupby('year').ha.sum().reindex(years).fillna(0).values; T['area_tc75_ha'] = ar[ar.tcb >= 3].groupby('year').ha.sum().reindex(years).fillna(0).values; T['area_tc25_ha'] = ar[ar.tcb == 0].groupby('year').ha.sum().reindex(years).fillna(0).values
CATS = {'颱風': 'Typhoon', '豪雨': 'Rainfall', '地震': 'Earthquake', '混合事件': 'Multiple typhoon and rainfall events'}; trig = CY.cat.map(CATS).fillna('Other or unspecified')
org = np.where(~hasg, 'no forest type nearby', np.where(L.func.values == 0, 'not mapped as forest stand', np.where(np.isin(L.func.values, [1, 2]), 'natural forest', np.where(L.func.values == 3, 'semi-natural forest', 'plantation'))))
szl = np.array(['<3 ha', '3-30 ha', '>=30 ha'])[szc]; ebl = np.array(['<500 m', '500-1500 m', '1500-2500 m', '>=2500 m'])[eb]
rows = []
for y, (i, l_) in lossy.items():
    cyy = CY[CY.year == y]
    for dim, lab in [('type', np.where(hasg, np.array(GS)[gi], 'none')[i]), ('origin', org[i]), ('size', szl[i]), ('elev', ebl[i]), ('trigger', trig[cyy.index].values), ('treecover2000', np.array(['0-25', '25-50', '50-75', '75-90', '90-100'])[ci[i]])]:
        s = pd.Series(l_).groupby(lab).sum() / 1e6; h = pd.Series(cyy.frac.values * 0.09).groupby(lab).sum()
        for kx in s.index: rows.append({'year': y, 'dim': dim, 'class': kx, 'loss_Mt': s[kx], 'area_ha': h[kx]})
B = pd.DataFrame(rows); B.to_csv(D + '/ledger_breakdown.csv', index=False)
for j, g in enumerate(GS): T['loss_' + g] = B[(B.dim == 'type') & (B['class'] == g)].set_index('year').loss_Mt.reindex(years).fillna(0).values
T.to_csv(D + '/ledger_annual.csv', index=False)
coh = pd.DataFrame({'first': first, 'loss': lossc, 'gain': gainc, 'deficit': (1 - f_end) * S_c * area, 'ha': area}).groupby('first').sum() / np.array([1e6, 1e6, 1e6, 1.0]); coh.to_csv(D + '/ledger_by_first_failure_year.csv')
L['loss_t'] = lossc; L['gx'] = (L.lon // 0.05) * 0.05 + 0.025; L['gy'] = (L.lat // 0.05) * 0.05 + 0.025; L['S_pre'] = S_c; L['area_ha'] = area
mp = L.groupby(['gx', 'gy']).agg(loss_t=('loss_t', 'sum'), ha=('area_ha', 'sum')).reset_index(); mp.to_csv(D + '/map_loss_0p05deg.csv', index=False)
blkl = L.groupby('blk').loss_t.sum().sort_values(ascending=False); nb_ = len(blkl)
par = {'x0': dict(zip(GS, x0.tolist())), 'x0_ci': dict(zip(GS, ci95(bx0).T.tolist())), 'k': dict(zip(GS, k.tolist())), 'k_ci': dict(zip(GS, ci95(bk).T.tolist())),
       'rho': {g: rho[j].tolist() for j, g in enumerate(GS)}, 'rho_ci': {g: ci95(brho)[:, j, :].T.tolist() for j, g in enumerate(GS)}, 'rho_ref': {g: rho_ref[j].tolist() for j, g in enumerate(GS)},
       'n_pre': {g: [int(((pre.tg == g) & (pre.tcb == c)).sum()) for c in range(5)] for g in GS}, 'n_post_age8_12_tc75': {g: int(((F.tg == g) & (F.age >= 8) & (F.age < 12)).sum()) for g in GS},
       'k_size6': {str(q): fe(F[(F.sz6 == q) & (F.tg != 'BB')], ['eb', 'tg']) for q in range(6)}, 'n_size6': {str(q): int(((F.sz6 == q) & (F.tg != 'BB')).sum()) for q in range(6)}, 'n_F': {g: int((F.tg == g).sum()) for g in GS},
       'n_post_age_lt1': {g: int(((post.tg == g) & (post.age < 1)).sum()) for g in GS}, 'k_specs': kspec, 'size_specific': {'x0': x0s.tolist(), 'k': ks.tolist()}, 'n_post': int(len(post)), 'n_pre_total': int(len(pre)), 'n_ref_local': int(len(Rf)),
       'href_mean': float(post.href.mean()), 'mean_Spre_tc90_BL': float(S_c[(gi == 0) & (ci == 4) & hasg].mean()), 'mean_Sbase_BL': float(S_base[(gi == 0) & hasg].mean())}
json.dump(par, open(D + '/parameters.json', 'w'), indent=1)
tot_t = B[B.dim == 'trigger'].groupby('class').loss_Mt.sum().to_dict()
summ = {'loss_total': float(loss_c.sum()), 'loss_total_ci': ci95(bl.sum(1)).tolist(), 'loss_2009': float(loss_c[5]), 'loss_2009_ci': ci95(bl[:, 5]).tolist(), 'gain_total': float(gain_c.sum()), 'gain_total_ci': ci95(bg.sum(1)).tolist(),
        'gain_2025': float(gain_c[-1]), 'gain_2025_ci': ci95(bg[:, -1]).tolist(), 'deficit_2025': float(def_c), 'deficit_ci': ci95(bd).tolist(), 'offset_pct': float(100 * gain_c.sum() / loss_c.sum()), 'offset_pct_ci': ci95(100 * bg.sum(1) / bl.sum(1)).tolist(),
        'loss_mean': float(loss_c.mean()), 'loss_median': float(np.median(loss_c)), 'top4_share': float(np.sort(loss_c)[-4:].sum() / loss_c.sum()), 'top_years': [int(y) for y in years[np.argsort(loss_c)[::-1][:5]]],
        'area_total_ha': float(area.sum()), 'area_tc75_ha': float(area[ci >= 3].sum()), 'area_tc25_ha': float(area[ci == 0].sum()), 'area_no_type_ha': float(area[~hasg].sum()), 'cells': int(len(L)), 'cells_2plus': int((L.nyears >= 2).sum()),
        'event_area_sum_ha': float((CY.frac * 0.09).sum()), 'nir_dist_2004_2024': float(sum(NIR_DIST.values()) / 1e3), 'loss_2004_2024': float(loss_c[:-1].sum()), 'nir_dist_2009': 2.535, 'nir_removal_mean': float(np.mean(list(NIR_REM.values())) / 1e3),
        'nir_area_2010_2024_ha': float(sum(NIR_AREA.values())), 'loss_by_type': B[B.dim == 'type'].groupby('class').loss_Mt.sum().to_dict(), 'loss_by_origin': B[B.dim == 'origin'].groupby('class').loss_Mt.sum().to_dict(), 'area_by_origin': B[B.dim == 'origin'].groupby('class').area_ha.sum().to_dict(),
        'loss_by_size': B[B.dim == 'size'].groupby('class').loss_Mt.sum().to_dict(), 'area_by_size': B[B.dim == 'size'].groupby('class').area_ha.sum().to_dict(), 'loss_by_elev': B[B.dim == 'elev'].groupby('class').loss_Mt.sum().to_dict(), 'loss_by_trigger': tot_t,
        'loss_by_tc': B[B.dim == 'treecover2000'].groupby('class').loss_Mt.sum().to_dict(), 'area_by_tc': B[B.dim == 'treecover2000'].groupby('class').area_ha.sum().to_dict(),
        'blocks_with_loss': int((blkl > 0).sum()), 'top10pct_blocks_share': float(blkl.iloc[:max(1, int(0.1 * (blkl > 0).sum()))].sum() / blkl.sum()), 'top20_blocks_share': float(blkl.iloc[:20].sum() / blkl.sum()),
        'variants': V, 'projection_deficit': proj, 'nboot': NB, 'cohort_2009': coh.loc[2009].to_dict(), 'mean_stock_lost_t_per_ha_2009': float(loss_c[5] * 1e6 / T.area_ha[5])}
json.dump(summ, open(D + '/ledger_summary.json', 'w'), indent=1, ensure_ascii=False)
np.savez_compressed(D + '/bootstrap.npz', loss=bl, gain=bg, deficit=bd, x0=bx0, k=bk, rho=brho)
L[['cell', 'a', 'first', 'last', 'nyears', 'tg', 'tcb', 'eb', 'sz', 'S_pre', 'loss_t', 'lon', 'lat']].to_parquet(D + '/cells_ledger.parquet')
pd.set_option('display.width', 250)
print('x0', np.round(x0, 3), 'k', np.round(k, 4), 'k CI', np.round(ci95(bk).T, 4).tolist()); print('rho'); print(np.round(rho, 3)); print('rho_ref'); print(np.round(rho_ref, 3))
print(T.round(3)[['year', 'area_ha', 'area_tc75_ha', 'loss', 'loss_lo', 'loss_hi', 'gain', 'gain_lo', 'gain_hi', 'nir_disturbance', 'loss_pct_removal']].to_string(index=False))
print(json.dumps({k_: v for k_, v in summ.items() if k_ not in ('variants',)}, ensure_ascii=False)[:3200]); print(pd.DataFrame(V).T.round(3).to_string()); print(json.dumps(kspec)[:900])
