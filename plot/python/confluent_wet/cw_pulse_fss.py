#!/usr/bin/env python3
"""Finite-size scan of the pulse response: recovery time and fluctuation.

Per case (called by submit_case.sh after the simulation):
    cw_pulse_fss.py RAW_CASE OUT_CASE --manifest MANIFEST
Summary (all sizes; the L256 reference runs are reduced from their raw dirs on first use):
    cw_pulse_fss.py RESULTS_ROOT OUT_SUMMARY --manifest MANIFEST --summary

Per case, response.csv is reduced to 1 tau_c bins whose edges sit on the pulse end
(t = 0), for every run of a pair alike, so pulse and control bins coincide exactly.

Two headline metrics, per (L, tau_m, init):
  T_half   half-recovery time: first t at which the seed-averaged pulse-minus-control
           <chi> response, smoothed over 5 tau_c, falls below half its value in the first
           2 tau_c after the pulse.  Model-free; equals tau ln2 for an exponential.
           (critical slowing down / recovery rate: van Nes & Scheffer 2007, Veraart et al.
           2012, Dai et al. 2012)
  sigma2   variance of <chi>(t) in the unperturbed controls over the stationary window
           t in [-1000, 1000] tau_c (the last 2000 tau_c of 3000 tau_c with feedback);
           S = L^2 sigma2 is the susceptibility, size-independent away from a critical point.
           (rising variance: Carpenter & Brock 2006, Scheffer et al. 2009)
Supporting: T_1/e, escape count (seed whose last-300 tau_c response still exceeds half
the initial response, same sign), integrated autocorrelation time tau_ac of the controls,
D_eff = sigma2 / tau_ac, Binder cumulant U = 1 - <d^4>/(3<d^2>^2), the across-seed
ensemble variance, and the same quantities for m. Errors: delete-one-seed jackknife.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np

TC = 674.3290333006435
COLUMNS = ('chi_mean', 'chi_std', 'm_mean', 'm_std', 'P_std', 'u_rms')
OBS = ('chi_mean', 'm_mean')
WINDOW = (-1000., 1000.)       # stationary window for control fluctuations, tau_c
LATE = 300.                    # escape test: last LATE tau_c of the recovery
SMOOTH = 5                     # running-mean width for crossing times, tau_c bins
ACF_MAX = 600                  # tau_c


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 23), b''):
            h.update(block)
    return h.hexdigest()


def parameters(root):
    data = json.loads((root / 'parameters.json').read_text())['data']
    return {k: v['value'] for k, v in data.items()}


def reduce_raw(root, pulse_end):
    """1 tau_c bins of response.csv, edges at pulse_end + k TC; unweighted bin means."""
    lines = '\n'.join(s for s in (root / 'response.csv').read_text().splitlines()
                      if not s.startswith('#'))
    d = np.atleast_1d(np.genfromtxt(io.StringIO(lines), delimiter=',', names=True))
    t = (d['step'] - pulse_end) / TC
    k = np.floor(t).astype(int)
    kmin, kmax = int(np.ceil(t.min())), int(np.floor(t.max()))   # fully covered bins only
    ok = (k >= kmin) & (k < kmax)
    idx = k[ok] - kmin
    n = np.bincount(idx, minlength=kmax - kmin)
    out = {'t_tc': np.arange(kmin, kmax) + .5, 'count': n}
    for c in COLUMNS:
        out[c] = np.bincount(idx, weights=d[c][ok], minlength=kmax - kmin) / np.maximum(n, 1)
    if np.any(n == 0):
        raise ValueError(f'Empty 1 tau_c bin in {root}')
    return out


def check_runtime(p, item):
    want = {'LX': item['L'], 'LY': item['L'], 'seed': item['seed'], 'chi_seed': item['chi_seed'],
            'chi0': int(item['initialization'])}
    for key, value in want.items():
        if float(p[key]) != float(value):
            raise ValueError(f'Runtime mismatch {key}: {p[key]} != {value}')
    if abs(float(p['tau_m']) / TC - item['tm_over_tc']) > 1e-4:
        raise ValueError('Runtime tau_m mismatch')
    steps = int(float(p['pmem_pulse_steps']))
    if (steps > 0) != (item['kind'] == 'pulse'):
        raise ValueError('Pulse/control identity mismatch')
    if int(float(p['pmem_pulse_start'])) + round(3 * TC) != item['pulse_end_steps']:
        raise ValueError('Pulse clock mismatch')


def reduce_one(raw, out, item, manifest_sha):
    p = parameters(raw)
    check_runtime(p, item)
    series = reduce_raw(raw, item['pulse_end_steps'])
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / 'fss_case.npz', **series)
    meta = dict(case=item['case'], raw=str(raw), manifest_sha256=manifest_sha,
                response_csv_sha256=sha(raw / 'response.csv'),
                analysis_sha256=sha(Path(__file__)), parameters=p)
    (out / 'fss_case.json').write_text(json.dumps(meta, indent=1, default=str) + '\n')


# ---------------------------------------------------------------- statistics
def jackknife(fn, seeds):
    """fn(list of seed arrays) -> scalar; returns (estimate, delete-one SE)."""
    full = fn(seeds)
    n = len(seeds)
    if n < 3 or not np.isfinite(full):
        return full, np.nan
    loo = np.array([fn(seeds[:i] + seeds[i + 1:]) for i in range(n)])
    if not np.isfinite(loo).all():
        return full, np.nan
    return full, float(np.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2)))


def crossing(t, r, level):
    """First t >= 0 at which the smoothed r/r0 falls below level."""
    post = t > 0
    t, r = t[post], r[post]
    r0 = r[:2].mean()
    if r0 == 0:
        return np.nan
    u = np.convolve(r / r0, np.ones(SMOOTH) / SMOOTH, mode='same')
    hit = np.nonzero(u[SMOOTH // 2:] < level)[0]
    return float(t[SMOOTH // 2:][hit[0]]) if len(hit) else np.nan


def acf_time(xs):
    """Integrated autocorrelation time (tau_c) of the seed-averaged normalized ACF,
    Sokal window: smallest W with W >= 5 tau_int(W)."""
    acs = []
    for x in xs:
        d = x - x.mean()
        f = np.fft.rfft(d, 2 * len(d))
        a = np.fft.irfft(f * np.conj(f))[:len(d)]
        acs.append(a / a[0])
    a = np.mean(acs, axis=0)[:ACF_MAX]
    tau = .5 + np.cumsum(a[1:])
    for w in range(1, len(tau)):
        if w >= 5 * tau[w - 1]:
            return float(tau[w - 1])
    return np.nan


def group_stats(L, pulses, controls):
    """pulses/controls: lists (same seed order) of reduced dicts."""
    res = {'n_seeds': len(pulses)}
    t = pulses[0]['t_tc']
    for c in pulses + controls:
        if not np.array_equal(c['t_tc'], t):
            raise ValueError('Seed/pair time grids differ')
    post = t > 0
    win = (t > WINDOW[0]) & (t < WINDOW[1])
    late = t > t[-1] - LATE
    for obs in OBS:
        r = [p[obs] - c[obs] for p, c in zip(pulses, controls)]
        r0 = float(np.mean([x[post][:2].mean() for x in r]))
        esc = [bool(np.sign(x[late].mean()) == np.sign(r0) and abs(x[late].mean()) > .5 * abs(r0))
               for x in r]
        o = {'r0': r0, 'escaped_seeds': int(sum(esc)),
             'late_over_r0_per_seed': [float(x[late].mean() / r0) for x in r]}
        for name, lev in (('T_half', .5), ('T_1e', np.exp(-1))):
            if any(esc):
                o[name], o[name + '_se'] = np.nan, np.nan
            else:
                o[name], o[name + '_se'] = jackknife(lambda s: crossing(t, np.mean(s, 0), lev), r)
        x = [c[obs][win] for c in controls]
        within = lambda s: float(np.mean([np.var(v, ddof=1) for v in s]))
        o['sigma2'], o['sigma2_se'] = jackknife(within, x)
        o['S'], o['S_se'] = L ** 2 * o['sigma2'], L ** 2 * o['sigma2_se']
        o['sigma2_ensemble'] = float(np.mean(np.var(np.array(x), axis=0, ddof=1)))
        o['seed_means'] = [float(v.mean()) for v in x]
        o['tau_ac'], o['tau_ac_se'] = jackknife(acf_time, x)
        o['D_eff'] = o['sigma2'] / o['tau_ac'] if o['tau_ac'] else np.nan

        def binder(s, pooled=False):
            d = np.concatenate([v - (np.mean(np.concatenate(s)) if pooled else v.mean()) for v in s])
            return float(1 - np.mean(d ** 4) / (3 * np.mean(d ** 2) ** 2))
        o['binder'], o['binder_se'] = jackknife(binder, x)
        o['binder_pooled'] = binder(x, True)
        o['skew'] = float(np.mean([np.mean((v - v.mean()) ** 3) / np.std(v) ** 3 for v in x]))
        half = len(x[0]) // 2
        o['drift_over_sigma'] = float(np.mean([(v[half:].mean() - v[:half].mean()) / v.std() for v in x]))
        res[obs] = o
    res['chi_std_spatial'] = float(np.mean([c['chi_std'][win].mean() for c in controls]))
    res['response_mean'] = {obs: np.mean([p[obs] - c[obs] for p, c in zip(pulses, controls)], 0)[post]
                            for obs in OBS}
    res['t_post'] = t[post]
    return res


# ---------------------------------------------------------------- summary
def load_reduced(path):
    with np.load(path) as d:
        return {k: d[k] for k in d.files}


def summarize(root, out, manifest_path):
    manifest = json.loads(manifest_path.read_text())
    msha = sha(manifest_path)
    out.mkdir(parents=True, exist_ok=True)
    items, missing = {}, []
    for item in manifest['cases']:
        f = root / item['case'] / 'fss_case.npz'
        if f.is_file():
            items[item['case']] = (item, load_reduced(f))
        else:
            missing.append(item['case'])
    for item in manifest['reference_L256']:
        target = out / 'reference' / item['case']
        if not (target / 'fss_case.npz').is_file():
            reduce_one(Path(item['raw']), target, item, msha)
        items[item['case']] = (item, load_reduced(target / 'fss_case.npz'))
    groups = {}
    for case, (item, data) in items.items():
        if item['kind'] != 'pulse':
            continue
        ctrl_case = case.replace('/pulse/', '/control/')
        if ctrl_case not in items:
            missing.append(ctrl_case)
            continue
        key = (item['L'], item['tm_over_tc'], str(item['initialization']))
        groups.setdefault(key, []).append((item['replicate'], data, items[ctrl_case][1]))
    rows, curves = [], {}
    for key in sorted(groups):
        members = sorted(groups[key], key=lambda m: m[0])
        L, g, init = key
        res = group_stats(L, [m[1] for m in members], [m[2] for m in members])
        name = f'L{L}_tm{g:g}_init{init}'.replace('.', 'p')
        curves[name + '_t'] = res.pop('t_post')
        for obs, v in res.pop('response_mean').items():
            curves[f'{name}_{obs}'] = v
        row = {'L': L, 'tm_over_tc': g, 'init': init, 'n_seeds': res['n_seeds'],
               'chi_std_spatial': res['chi_std_spatial']}
        for obs in OBS:
            for k, v in res[obs].items():
                row[f'{obs}.{k}'] = v
        rows.append(row)
        print(json.dumps({'group': name, 'seeds': res['n_seeds'],
                          'T_half': res['chi_mean']['T_half'], 'S': res['chi_mean']['S'],
                          'escaped': res['chi_mean']['escaped_seeds']}), flush=True)
    np.savez_compressed(out / 'fss_response_curves.npz', **curves)
    result = dict(manifest=str(manifest_path), manifest_sha256=msha,
                  analysis_sha256=sha(Path(__file__)), missing=sorted(set(missing)),
                  window_tc=WINDOW, late_tc=LATE, smooth_tc=SMOOTH, groups=rows,
                  definitions=__doc__)
    (out / 'fss_summary.json').write_text(json.dumps(result, indent=1, default=float) + '\n')
    flat = [{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()} for r in rows]
    with (out / 'fss_groups.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, list(flat[0]))
        w.writeheader()
        w.writerows(flat)
    figures(out, rows)
    print(json.dumps({'groups': len(rows), 'missing': len(set(missing))}), flush=True)
    return not missing


def figures(out, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {128: '#D68910', 256: '#2471A3', 512: '#7D3C98'}
    panels = (('chi_mean.T_half', 'half-recovery time T_half / tau_c', False),
              ('chi_mean.S', 'S = L^2 Var<chi>', True),
              ('chi_mean.tau_ac', 'autocorrelation time tau_ac / tau_c', False),
              ('chi_mean.binder', 'Binder cumulant U (within run)', False))
    for init in ('0', '1'):
        fig, axes = plt.subplots(1, 4, figsize=(17, 3.8))
        for ax, (key, label, logy) in zip(axes, panels):
            for L in sorted({r['L'] for r in rows}):
                sel = sorted((r for r in rows if r['L'] == L and r['init'] == init),
                             key=lambda r: r['tm_over_tc'])
                x = [r['tm_over_tc'] for r in sel]
                y = [r[key] for r in sel]
                e = [r.get(key + '_se', np.nan) for r in sel]
                ax.errorbar(x, y, e, fmt='o-', color=colors.get(L, 'k'), mfc='none', ms=4,
                            capsize=2, lw=1, label=f'L={L}')
                esc = [r for r in sel if r['chi_mean.escaped_seeds'] > 0]
                if key.endswith('T_half') and esc:
                    ax.plot([r['tm_over_tc'] for r in esc], [5] * len(esc), 'x', color=colors.get(L, 'k'))
            ax.set_xlabel('tau_m / tau_c')
            ax.set_title(label, fontsize=10)
            if logy:
                ax.set_yscale('log')
        axes[0].legend(fontsize=8)
        axes[0].text(.02, .02, 'x: >=1 seed escaped (no T_half)', transform=axes[0].transAxes, fontsize=7)
        fig.suptitle(f'Finite-size scan, init {init}: recovery and fluctuation of <chi>')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(out / f'fss_init{init}.{ext}', dpi=150)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root', type=Path)
    ap.add_argument('out', type=Path)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--summary', action='store_true')
    ap.add_argument('--allow-missing', action='store_true',
                    help='summary of whatever is complete (e.g. the L256 reference alone)')
    a = ap.parse_args()
    manifest_path = a.manifest.resolve()
    if a.summary:
        if not summarize(a.root, a.out, manifest_path) and not a.allow_missing:
            raise SystemExit(2)
        return
    manifest = json.loads(manifest_path.read_text())
    raw = a.root.absolute()
    match = [c for c in manifest['cases'] if raw.as_posix().endswith('/' + c['case'])]
    if len(match) != 1:
        raise ValueError(f'Cannot identify the manifest case of {raw}')
    item = match[0]
    card = manifest_path.parent / item['case'] / 'run.dat'
    if sha(card) != item['run_dat_sha256']:
        raise ValueError('Tracked run.dat differs from manifest')
    reduce_one(raw, a.out, item, sha(manifest_path))


if __name__ == '__main__':
    main()
