#!/usr/bin/env python3
"""Direct p1/p2 comparison. Per-case postprocessing and a dependent scan summary."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[3]
MANIFEST = REPO / 'cases/20260915/cw_pressure_compare/manifest.json'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(part)
    return h.hexdigest()


def read(path):
    data = json.loads(path.read_text())['data']
    return {k: v['value'] for k, v in data.items()}


def analyse(root, out):
    manifest = json.loads(MANIFEST.read_text())
    row = next(r for r in manifest['cases'] if root.as_posix().endswith('/'+r['case']))
    par = read(root/'parameters.json')
    # Check every serialized counterpart of the input; tolerate writer's six-digit precision.
    for key, expected in row['expected_parameters'].items():
        actual_key = key.replace('-', '_')
        if actual_key not in par:
            continue
        actual = par[actual_key]
        try:
            ok = np.isclose(float(actual), float(expected), rtol=6e-6, atol=1e-12)
        except (TypeError, ValueError):
            ok = str(actual) == expected
        if not ok:
            raise ValueError(f'Runtime mismatch {key}: {actual} != {expected}')
    for required in ('LX', 'LY', 'nsteps', 'ninfo', 'open_loop', 'zeta_open', 'frame_light'):
        if required not in par:
            raise ValueError(f'Missing required runtime parameter: {required}')
    p1s, p2s, lbs, trace, hashes = [], [], [], [], {}
    max_error = 0.
    for t in range(0, manifest['nsteps']+1, manifest['ninfo']):
        path = root/f'frame{t}.json'
        data = read(path)
        p1 = np.asarray(data['pressure'], dtype=float).ravel()
        p2 = -np.asarray(data['sigma_bulk'], dtype=float).ravel()
        lb = np.asarray(data['pressure_lb'], dtype=float).ravel()
        ux = np.asarray(data['ux_mat'], dtype=float).ravel()
        uy = np.asarray(data['uy_mat'], dtype=float).ravel()
        for arr in (p1, p2, lb, ux, uy):
            if arr.size != 256**2 or not np.isfinite(arr).all():
                raise ValueError(f'Invalid native field at {t}')
        err = float(np.max(np.abs(p1-p2-lb)))
        max_error = max(err, max_error)
        if not np.allclose(p1, lb+p2, rtol=2e-5, atol=2e-7):
            raise ValueError(f'Pressure decomposition failed at {t}: {err}')
        if np.max(p2) > 0:
            raise ValueError('p2 must be nonpositive')
        trace.append([t, np.sqrt(np.mean(ux*ux+uy*uy)), p1.mean(), p1.std(),
                      p2.mean(), p2.std(), lb.mean(), lb.std()])
        hashes[path.name] = sha(path)
        if t >= manifest['steady_start']:
            p1s.append(p1); p2s.append(p2); lbs.append(lb)
    p1, p2, lb = map(np.concatenate, (p1s, p2s, lbs))
    trace = np.asarray(trace)
    window = trace[trace[:, 0] >= manifest['steady_start']]
    first, second = np.array_split(window, 2)
    drift = {name: float(abs(first[:, col].mean()-second[:, col].mean()) /
                         max(abs(window[:, col].mean()), 1e-30))
             for name, col in [('u_rms', 1), ('p1_sigma', 3), ('p2_sigma', 5)]}
    for label, col, arr in [('p1', 2, p1), ('p2', 4, p2)]:
        drift[label+'_mean_in_sigma'] = float(abs(first[:, col].mean()-second[:, col].mean()) /
                                                max(arr.std(), 1e-30))
    cov = float(np.mean((lb-lb.mean())*(p2-p2.mean())))
    stats = {f'{label}_{metric}': float(fun(arr)) for label, arr in [('p1', p1), ('p2', p2), ('pLB', lb)]
             for metric, fun in [('mean', np.mean), ('sigma', np.std)]}
    result = dict(case=row['case'], activity=row['activity'], seed=row['seed'],
                  replicate=row['replicate'], zeta_open=row['zeta_open'], **stats,
                  covariance_pLB_p2=cov, variance_identity_residual=float(p1.var()-lb.var()-p2.var()-2*cov),
                  stationary=bool(max(drift.values()) <= .15), drift=drift,
                  steady_frames=len(window), sites_per_frame=256**2,
                  steady_steps=[int(window[0, 0]), int(window[-1, 0])],
                  max_decomposition_error=max_error, frame_sha256=hashes,
                  parameters_sha256=sha(root/'parameters.json'),
                  manifest_sha256=sha(MANIFEST), analysis_sha256=sha(Path(__file__)),
                  run_dat_sha256=row['run_dat_sha256'])
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out/'pressure_samples.npz', p1=p1, p2=p2)
    np.savetxt(out/'time_series.csv', trace, delimiter=',',
               header='step,u_rms,p1_mean,p1_sigma,p2_mean,p2_sigma,pLB_mean,pLB_sigma', comments='')
    (out/'pressure_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'frame_sha256'}))


def save(fig, out, name):
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'{name}.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


def summarize(root, out):
    manifest = json.loads(MANIFEST.read_text())
    out.mkdir(parents=True, exist_ok=True)
    rows, missing = [], []
    for case in manifest['cases']:
        path = root/case['case']/'pressure_result.json'
        if not path.exists():
            missing.append(case['case']); continue
        row = json.loads(path.read_text())
        if row['manifest_sha256'] != sha(MANIFEST) or row['run_dat_sha256'] != case['run_dat_sha256']:
            raise ValueError(f'Provenance mismatch: {path}')
        rows.append(row)
    report = dict(expected=len(manifest['cases']), available=len(rows), missing=missing,
                  nonstationary=[r['case'] for r in rows if not r['stationary']], cases=rows,
                  sigma_definition=manifest['sigma_definition'])
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    if missing:
        raise RuntimeError(f'Missing cases; no incomplete trend generated: {missing}')
    columns = ['case', 'activity', 'zeta_open', 'seed', 'stationary', 'p1_mean', 'p2_mean',
               'pLB_mean', 'p1_sigma', 'p2_sigma', 'pLB_sigma', 'covariance_pLB_p2',
               'variance_identity_residual']
    with (out/'pressure_vs_activity.csv').open('w') as f:
        w = csv.DictWriter(f, columns, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    activities = sorted({r['activity'] for r in rows})
    colors = {'p1': '#2166ac', 'p2': '#b2182b'}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    for label in colors:
        for ax, metric in zip(axes, ['mean', 'sigma']):
            values = [np.array([r[f'{label}_{metric}'] for r in rows if r['activity'] == a]) for a in activities]
            ax.errorbar(activities, [v.mean() for v in values], yerr=[v.std(ddof=1) for v in values],
                        marker='o', color=colors[label], capsize=3, label=label)
            for r in rows:
                ax.scatter(r['activity'], r[f'{label}_{metric}'], s=20, alpha=.5,
                           marker='o' if r['stationary'] else 'x', color=colors[label])
            ax.set(xlabel=r'$a=\zeta_{\mathrm{eff}}/\zeta$', ylabel=r'$\langle p\rangle$' if metric=='mean' else r'$\sigma_p$')
    axes[0].legend(frameon=False)
    save(fig, out, 'pressure_moments')
    fig, axes = plt.subplots(2, 3, figsize=(11, 6))
    separate = {label: plt.subplots(2, 3, figsize=(11, 6)) for label in colors}
    for a, ax, i in zip(activities, axes.flat, range(len(activities))):
        group = [r for r in rows if r['activity'] == a]
        samples = {label: [] for label in colors}
        for r in group:
            with np.load(root/r['case']/'pressure_samples.npz') as data:
                for label in colors:
                    samples[label].append(data[label])
        hist_export = {}
        for label in colors:
            data = np.concatenate(samples[label])
            density, edges = np.histogram(data, bins=240, density=True)
            hist_export[label+'_edges'] = edges; hist_export[label+'_density'] = density
            centers = (edges[1:]+edges[:-1])/2
            ax.plot(centers, density, color=colors[label], label=label)
            ax2 = separate[label][1].flat[i]
            ax2.plot(centers, density, color=colors[label])
            ax2.set(title=f'a={a:g}', xlabel=label, ylabel='Probability density')
            del data
        flagged = any(not r['stationary'] for r in group)
        ax.set(title=f'a={a:g}'+(' (drift flagged)' if flagged else ''), xlabel='p', ylabel='Probability density')
        np.savez_compressed(out/f'pdf_a{a:g}.npz', **hist_export)
        del samples
    axes.flat[0].legend(frameon=False)
    save(fig, out, 'pressure_distributions_overlay')
    for label, (fig, _) in separate.items():
        save(fig, out, f'pressure_distributions_{label}')
    fig, axes = plt.subplots(3, 2, figsize=(10, 8), sharex=True)
    for row in rows:
        data = np.loadtxt(root/row['case']/'time_series.csv', delimiter=',', skiprows=1)
        for ax, col, name in zip(axes.flat, [1, 2, 3, 4, 5, 6],
                                  ['u_rms', 'mean p1', 'sigma p1', 'mean p2', 'sigma p2', 'mean pLB']):
            ax.plot(data[:,0]/manifest['tau_c'], data[:,col], lw=.65, alpha=.6)
            ax.set(ylabel=name, xlabel=r'$t/\tau_c$')
            ax.axvline(manifest['steady_start']/manifest['tau_c'], color='k', ls=':', lw=.7) if row is rows[0] else None
    save(fig, out, 'stationarity_traces')
    (out/'README.md').write_text('# Pressure definition comparison\n\n'
        'p1 = pressure = pLB - sigma_bulk; p2 = -sigma_bulk. All three fields are read directly '
        'from full native JSON, with the same resting-state gauge; no centering of PDF inputs.\n\n'
        'Open-loop L256, 6 activities, 3 paired seeds each. Window: steps 134800..404400; '
        '41 equally weighted snapshots per seed. PDFs pool equal numbers of sites, snapshots and seeds; '
        '240 bins spanning the entire observed range for each definition and activity, no clipping or KDE. '
        'Separate p1/p2 plots use their own axes to resolve the narrow p2 distribution.\n\n'
        'Trend: population SD over space/time per seed; mean and sample SD across 3 seeds '
        '(error bars are seed spread, not confidence intervals). Individual seeds are shown; x means '
        'stationarity drift flagged. Two-half drift gate is a diagnostic, not proof of stationarity. '
        'See summary.json for flags and covariance/variance decomposition.\n\n'
        'p2 is nonpositive for CC>0. This study compares observables under prescribed activity; '
        'it does not substitute p2 into the closed-loop memory response or recalibrate its threshold.\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'cases'}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path); p.add_argument('out', type=Path)
    p.add_argument('--summary', action='store_true')
    args = p.parse_args()
    (summarize if args.summary else analyse)(args.root, args.out)
