#!/usr/bin/env python3
"""+1/2 motion relative to the moving active | passive interface, from dense npz runs (cw_frames_to_npz.py).

  cw_defect_interface_relative.py case <run.npz> <out_dir>     per-case statistics -> <out_dir>/<case>.json
  cw_defect_interface_relative.py plot <out_dir>               figure from every <case>.json in out_dir

Per frame: +1/2 cores by plaquette winding, head h = -div Q/|div Q|; interface = chi smoothed (Gaussian
sigma 1, periodic) at 1/2; signed distance s(x) (+ inside active, - inside passive; periodic EDT on a
3 x 3 tiling). Tracking: mutual nearest neighbours between consecutive frames within 2.5 grid. For a
+1/2 tracked from frame i to i+2:
  ds/dt = [s_{i+2}(x_{i+2}) - s_i(x_i)]/dt       motion relative to the interface (+ = deeper into active)
        = D + G,  D = [s_{i+2}(x_{i+2}) - s_{i+2}(x_i)]/dt   the defect's own displacement across the interface
                  G = [s_{i+2}(x_i) - s_i(x_i)]/dt           interface motion at the old position (+ = active grows)
Classes at frame i: HA head active / tail passive, HP head passive / tail active (mean raw chi along
core +- s h, s = 3..7, threshold .5), within |s_i(x_i)| <= 6 grid. Units grid per reference tau_c.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates
from scipy.spatial import cKDTree

TC = 674.3290333006435
BAND = 6.0
S_RAY = np.arange(3.0, 7.0 + 1e-9, 0.5)


def grad(f):
    return (np.roll(f, -1, 0) - np.roll(f, 1, 0))/2, (np.roll(f, -1, 1) - np.roll(f, 1, 1))/2


def plus_defects(qxx, qxy):
    th = .5*np.arctan2(qxy, qxx)
    d = lambda a, b: (b - a + np.pi/2) % np.pi - np.pi/2
    t10 = np.roll(th, -1, 0); t11 = np.roll(t10, -1, 1); t01 = np.roll(th, -1, 1)
    w = (d(th, t10) + d(t10, t11) + d(t11, t01) + d(t01, th))/(2*np.pi)
    return np.argwhere(w > .25) + .5


def at(f, xy):
    xy = np.asarray(xy); s = xy.shape[:-1]
    return map_coordinates(f, xy.reshape(-1, 2).T, order=1, mode='wrap').reshape(s)


def link(a, b, L, r=2.5):
    out = -np.ones(len(a), int)
    if not len(a) or not len(b):
        return out
    da, ia = cKDTree(b, boxsize=L).query(a)
    _, ib = cKDTree(a, boxsize=L).query(b)
    ok = (da < r) & (ib[ia] == np.arange(len(a)))
    out[ok] = ia[ok]
    return out


def signed_distance(chi):
    m = gaussian_filter(chi, 1, mode='wrap') > .5
    if m.all() or not m.any():
        return np.full(chi.shape, np.inf if m.all() else -np.inf)
    n = chi.shape[0]; t = np.tile(m, (3, 3))
    d_in = distance_transform_edt(t)[n:2*n, n:2*n]
    d_out = distance_transform_edt(~t)[n:2*n, n:2*n]
    return np.where(m, d_in - .5, -(d_out - .5))


def frame(z, i):
    qxx, qxy = z['QQxx'][i].astype(float), z['QQyx'][i].astype(float)
    chi = z['chi'][i].astype(float)
    plus = plus_defects(qxx, qxy)
    a, b = grad(qxx); c, e = grad(qxy)
    px, py = at(a + e, plus), at(c - b, plus)
    h = -np.c_[px, py]/np.maximum(np.hypot(px, py), 1e-300)[:, None]
    return plus, h, chi


def case(npz, out):
    z = np.load(npz); T = len(z['t']); L = z['chi'].shape[1]
    dt = 2*float(z['t'][1] - z['t'][0])/TC
    F = [frame(z, i) for i in range(T)]
    S = [signed_distance(f[2]) for f in F]
    links = [link(F[i][0], F[i + 1][0], L) for i in range(T - 1)]
    rows = []
    for i in range(T - 2):
        j1 = links[i]; ok = j1 >= 0
        j2 = np.where(ok, links[i + 1][np.maximum(j1, 0)], -1); ok &= j2 >= 0
        if not ok.any():
            continue
        plus, h, chi = F[i]
        x0, h0 = plus[ok], h[ok]; x2 = F[i + 2][0][j2[ok]]
        s0, s2n, s2o = at(S[i], x0), at(S[i + 2], x2), at(S[i + 2], x0)
        hh = at(chi, x0[:, None, :] + S_RAY[None, :, None]*h0[:, None, :]).mean(1)
        tt = at(chi, x0[:, None, :] - S_RAY[None, :, None]*h0[:, None, :]).mean(1)
        v = np.sum(((x2 - x0 + L/2) % L - L/2)*h0, 1)/dt
        rows.append(np.c_[s0, (s2n - s0)/dt, (s2n - s2o)/dt, (s2o - s0)/dt, hh, tt, v])
    A = np.concatenate(rows); A = A[np.isfinite(A).all(1)]
    s0, rel, D, G, hh, tt, v = A.T
    near = np.abs(s0) <= BAND
    ha, hp = (hh > .5) & (tt <= .5), (hh <= .5) & (tt > .5)
    stat = lambda x, sel: {'mean': float(x[sel].mean()), 'se': float(x[sel].std()/np.sqrt(max(sel.sum(), 1)))}
    res = {k: {'n': int(sel.sum())} | {q: stat(x, sel) for q, x in (('ds_dt', rel), ('D', D), ('G', G), ('v_head', v))}
           for k, sel in (('HA', near & ha), ('HP', near & hp), ('all_near', near)) if sel.sum()}
    edges = np.arange(-BAND, BAND + 1e-9, 1.0); b = np.digitize(s0, edges) - 1
    prof = {}
    for k, cl in (('HA', ha), ('HP', hp)):
        keep = [j for j in range(len(edges) - 1) if (cl & (b == j)).sum() >= 100]
        prof[k] = {'s': [float(s0[cl & (b == j)].mean()) for j in keep], 'n': [int((cl & (b == j)).sum()) for j in keep]}
        for q, x in (('ds_dt', rel), ('D', D), ('G', G)):
            prof[k][q] = [float(x[cl & (b == j)].mean()) for j in keep]
            prof[k][q + '_se'] = [float(x[cl & (b == j)].std()/np.sqrt((cl & (b == j)).sum())) for j in keep]
    name = Path(npz).stem
    out.mkdir(parents=True, exist_ok=True)
    (out/f'{name}.json').write_text(json.dumps({'case': name, 'dt_tc': dt, 'samples': int(len(s0)), 'summary': res,
                                                'profiles': prof}, indent=1))
    print(name, {k: (r['n'], round(r['ds_dt']['mean'], 3), round(r['D']['mean'], 3), round(r['G']['mean'], 3)) for k, r in res.items()})


def plot(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    J = [json.loads(f.read_text()) for f in sorted(out.glob('hn*.json'))]
    J = [j for j in J if 'HA' in j['summary'] and 'HP' in j['summary']]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False})
    fig, ax = plt.subplots(2, len(J), figsize=(5.2*len(J), 8), squeeze=False)
    for j, r in enumerate(J):
        x = np.arange(3); w = .38
        for k, (cl, col, lab) in enumerate((('HA', '#c8382f', 'head active | tail passive'), ('HP', '#1f5fa8', 'head passive | tail active'))):
            s = r['summary'][cl]
            ax[0, j].bar(x + (k - .5)*w, [s[q]['mean'] for q in ('D', 'G', 'ds_dt')], w,
                         yerr=[s[q]['se'] for q in ('D', 'G', 'ds_dt')], color=col, capsize=2, label=f'{lab} (n {s["n"]})')
            p = r['profiles'][cl]
            ax[1, j].errorbar(p['s'], p['ds_dt'], yerr=p['ds_dt_se'], fmt='-o', ms=4, color=col, capsize=2, label=f'{cl}: relative ds/dt')
            ax[1, j].plot(p['s'], p['D'], '--', color=col, lw=1, label=f'{cl}: defect motion D')
            ax[1, j].plot(p['s'], p['G'], ':', color=col, lw=1.3, label=f'{cl}: interface motion G')
        ax[0, j].axhline(0, color='0.5', lw=.8)
        ax[0, j].set_xticks(x); ax[0, j].set_xticklabels(['defect motion D', 'interface motion G\n(+ active grows)', 'relative ds/dt = D + G'])
        ax[0, j].set_ylabel('normal velocity (grid/$\\tau_c$), + = into active'); ax[0, j].set_title(r['case'], fontsize=9)
        ax[0, j].legend(fontsize=8)
        ax[1, j].axhline(0, color='0.5', lw=.8); ax[1, j].axvline(0, color='0.5', lw=.8)
        ax[1, j].set_xlabel('signed distance to interface s (grid, + inside active)'); ax[1, j].set_ylabel('grid/$\\tau_c$')
        ax[1, j].legend(fontsize=7.5)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'interface_relative.{ext}', dpi=200)
    print('figure', out/'interface_relative.png', len(J), 'cases')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode', choices=('case', 'plot'))
    ap.add_argument('paths', nargs='+', type=Path)
    a = ap.parse_args()
    if a.mode == 'case':
        case(a.paths[0], a.paths[1])
    else:
        plot(a.paths[0])
