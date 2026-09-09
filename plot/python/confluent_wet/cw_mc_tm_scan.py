#!/usr/bin/env python3
"""Terminal-state diagnostics and empirical phase map for cw_mc_tm_phase.

Per run: INPUT_RAW OUTPUT_RESULTS (the submit_case.sh PLOT_SCRIPT contract).
Whole campaign: INPUT_RESULTS OUTPUT_SUMMARY --summary --manifest manifest.json.
Only exact full-resolution domain averages from video_meta.csv are used; no
mean is estimated from uint8 pixels. Phase labels describe finite-time evidence.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re

import numpy as np

TC = 674.3290333006435
DRIFT = 0.03
BLOCK_RANGE = 0.06


def parameters(root):
    d = json.loads((root/'parameters.json').read_text())['data']
    return {k: v['value'] for k, v in d.items()}


def block_means(t, x, start, end, n):
    edges = np.linspace(start, end, n+1)
    return np.array([np.mean(x[(t >= a) & ((t < b) if i < n-1 else (t <= b))])
                     for i, (a, b) in enumerate(zip(edges[:-1], edges[1:]))])


def reduce_case(root, out):
    p = parameters(root)
    text = '\n'.join(s for s in (root/'video_meta.csv').read_text().splitlines()
                     if not s.startswith('#'))
    d = np.atleast_1d(np.genfromtxt(io.StringIO(text), delimiter=',', names=True))
    if len(d) < 20 or any(not np.isfinite(d[k]).all() for k in d.dtype.names):
        raise ValueError('Missing or non-finite domain-average record')
    if np.any(np.diff(d['t']) != int(p['nvideo'])):
        raise ValueError('Irregular or duplicate stream timestamps')
    if np.min(d['chi_mean']) < -1e-6 or np.max(d['chi_mean']) > 1+1e-6:
        raise ValueError('Mean chi outside [0,1]')
    nsteps = int(p['nsteps']); freeze = int(p.get('chi_freeze_steps', 0))
    # parameters.json uses the C++ stream's limited decimal precision. Keep the
    # requested grid coordinate from the case path, and check the serialized value.
    tm_token = re.search(r'_tm([0-9]+(?:p[0-9]+)?)_chi', root.name)
    if tm_token is None: raise ValueError('Missing memory-time coordinate in case name')
    tm = float(tm_token.group(1).replace('p', '.'))
    tm_serialized = float(p['tau_m'])/TC
    if not np.isclose(tm_serialized, tm, rtol=5e-6):
        raise ValueError('Case path and serialized tau_m disagree')
    t = (d['t']-freeze)/TC
    end = (nsteps-freeze)/TC
    window = max(500., 20*tm)
    if end < 2*window:
        raise ValueError('Observation too short for the campaign terminal windows')
    complete = (d['t'][-1] == (nsteps//int(p['nvideo']))*int(p['nvideo'])
                and (root/f'frame{nsteps}.json').is_file())
    tail = (t >= end-window) & (t <= end)
    if tail.sum() < 100:
        raise ValueError('Insufficient terminal samples')
    x = d['chi_mean']
    blocks = block_means(t, x, end-window, end, 4)
    halves = block_means(t, x, end-window, end, 2)
    # Slow excursions between both extreme phases are retained as a separate flag.
    nb = max(2, int(window/max(10., tm)))
    residence = block_means(t, x, end-window, end, nb)
    both_phases = bool(np.any(residence < 0.2) and np.any(residence > 0.8))
    drift = float(halves[1]-halves[0])
    flags = []
    if not complete: flags.append('incomplete_output')
    if abs(drift) > DRIFT: flags.append('terminal_drift')
    if np.ptp(blocks) > BLOCK_RANGE: flags.append('terminal_block_variation')
    if both_phases: flags.append('visits_both_phases')
    prep_window = min(freeze/TC/2, max(30., 3*tm))
    prep_m = block_means(t, d['m_mean'], -prep_window, 0., 2)
    prep_u = block_means(t, d['u_rms'], -prep_window, 0., 2)
    if (not np.isfinite(prep_m).all() or not np.isfinite(prep_u).all()
            or abs(prep_m[1]-prep_m[0]) > 0.01
            or abs(prep_u[1]-prep_u[0])/max(float(np.mean(prep_u)), 1e-12) > 0.15):
        flags.append('preparation_drift')
    if min(prep_u) <= 1e-8: flags.append('flow_not_developed')
    suffix = re.search(r'_rep([12])$', root.name)
    if suffix is None: raise ValueError('Missing replicate identifier in case name')
    row = {'case': f'L{int(p["LX"])}/{root.name}', 'L': int(p['LX']),
           'mc': float(p['mc']), 'tm_over_tc': tm, 'chi0': int(p['chi0']),
           'tm_over_tc_serialized': tm_serialized,
           'replicate': int(suffix.group(1)), 'seed': int(p['seed']),
           'tau_c_reference': TC, 'pmem': float(p['pmem']),
           'tau_chi': float(p['tau_chi']), 'nsteps': nsteps,
           'preparation_steps': freeze, 'observation_tc': end,
           'window_tc': window, 'window_tau_m': window/tm,
           'tail_mean': float(np.mean(x[tail])), 'tail_std_time': float(np.std(x[tail])),
           'tail_last': float(x[-1]), 'tail_half_drift': drift,
           'tail_four_blocks': blocks.tolist(), 'tail_block_range': float(np.ptp(blocks)),
           'coarse_residence_means': residence.tolist(),
           'tail_active_occupancy': float(np.mean(x[tail] < 0.2)),
           'tail_passive_occupancy': float(np.mean(x[tail] > 0.8)),
           'tail_m_mean': float(np.mean(d['m_mean'][tail])),
           'tail_pressure_std_mean': float(np.mean(d['P_std'][tail])),
           'tail_u_rms_mean': float(np.mean(d['u_rms'][tail])),
           'prep_m_halves': prep_m.tolist(), 'prep_u_rms_halves': prep_u.tolist(),
           'complete': bool(complete), 'diagnostic_flags': flags,
           'settled_by_diagnostics': not flags,
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'interpretation': 'Window-conditioned diagnostics; not a stationarity proof or iid confidence interval.'}
    if not np.isfinite(blocks).all(): raise ValueError('Empty terminal block')
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out/'series.npz', observation_time_tc=t,
                        **{k: d[k] for k in d.dtype.names})
    (out/'phase_case.json').write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'case': row['case'], 'tail_mean': row['tail_mean'], 'flags': flags}))


def classify(rows):
    if len(rows) != 4:
        return 'missing', None
    if any(not r['settled_by_diagnostics'] for r in rows):
        return 'unresolved', None
    arms = [[r['tail_mean'] for r in sorted(rows, key=lambda r: r['replicate'])
             if r['chi0'] == chi] for chi in (0, 1)]
    if any(len(a) != 2 for a in arms): return 'missing', None
    a, b = map(np.asarray, arms)
    if max(np.ptp(a), np.ptp(b)) > 0.08: return 'seed_sensitive', None
    if np.min(b-a) > 0.20 and np.min(b)-np.max(a) > 0.10:
        return 'initialization_dependent', None
    if np.ptp(np.r_[a, b]) <= 0.05:
        return 'initialization_independent', float(np.mean(np.r_[a, b]))
    return 'unresolved', None


def aggregate(root, out, manifest_path):
    manifest = json.loads(manifest_path.read_text())
    rows = []; missing = []; invalid = []
    for item in manifest['cases']:
        path = root/item['case']/'phase_case.json'
        if not path.exists():
            missing.append(item['case']); continue
        try:
            row = json.loads(path.read_text())
            for key in ('case', 'L', 'mc', 'chi0', 'replicate', 'seed', 'nsteps', 'preparation_steps'):
                if row[key] != item[key]: raise ValueError(f'Wrong {key}')
            if not np.isclose(row['tm_over_tc'], item['tm_over_tc'], rtol=1e-10):
                raise ValueError('Wrong tau_m')
            for key in ('pmem', 'tau_chi'):
                if row[key] != manifest['fixed'][key]: raise ValueError(f'Wrong {key}')
            rows.append(row)
        except (KeyError, ValueError) as exc:
            invalid.append({'case': item['case'], 'error': str(exc)})
    groups = []
    for Lname, grid in manifest['grids'].items():
        L = int(Lname[1:])
        for tm in grid['tm']:
            for mc in grid['mc']:
                rr = [r for r in rows if r['L'] == L and r['mc'] == mc
                      and np.isclose(r['tm_over_tc'], tm, rtol=1e-10)]
                state, value = classify(rr)
                arms = [[r['tail_mean'] for r in rr if r['chi0'] == a] for a in (0, 1)]
                groups.append({'L': L, 'tm_over_tc': tm, 'mc': mc,
                               'classification': state, 'unique_chi': value,
                               'chi0_means': arms[0], 'chi1_means': arms[1],
                               'case_count': len(rr),
                               'flags': {r['case']: r['diagnostic_flags'] for r in rr
                                         if r['diagnostic_flags']}})
    summary = {'expected_cases': len(manifest['cases']), 'available_cases': len(rows),
               'missing': missing, 'invalid': invalid, 'groups': groups,
               'criteria': manifest['classification'], 'cases': rows,
               'interpretation': manifest['interpretation']}
    out.mkdir(parents=True, exist_ok=True)
    (out/'phase_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    phase_figures(groups, manifest, out)
    counts = {k: sum(g['classification'] == k for g in groups)
              for k in ('initialization_independent', 'initialization_dependent',
                        'unresolved', 'seed_sensitive', 'missing')}
    (out/'README.md').write_text(
        '# Empirical phase map\n\n'
        f'Available {len(rows)}/{len(manifest["cases"])} cases. Classification counts: {counts}.\n\n'
        'Each point has two independent flow seeds and chi=0/1 starts. Displayed values '
        'are temporal means of exact full-field means over max(500 tau_c,20 tau_m). '
        'Preparation is excluded. Solid fill is assigned only when terminal and preparation '
        'diagnostics pass. Purple marks reproducible initialization dependence over the '
        'finite observation horizon; gray marks unresolved or seed-sensitive outcomes; '
        'x marks missing data. The two starting-arm means are shown separately. '
        'Point markers show actual sampled coordinates; cell colours do not interpolate a phase boundary.\n\n'
        'No iid significance test is applied to correlated time samples, and two seeds do '
        'not establish a thermodynamic or infinite-time phase boundary. Reference tau_c and '
        'physical coefficients are held fixed between sizes. Complete numeric diagnostics, '
        'flags and missing/invalid cases are in phase_summary.json.\n')
    print(json.dumps({'available': len(rows), 'expected': len(manifest['cases']),
                      'invalid': len(invalid), 'classifications': counts}))


def phase_figures(groups, manifest, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, Normalize
    from matplotlib.patches import Patch
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'pdf.fonttype': 42, 'axes.spines.top': False,
                         'axes.spines.right': False})
    for Lname, grid in manifest['grids'].items():
        L = int(Lname[1:]); tms = np.array(grid['tm']); mcs = np.array(grid['mc'])
        tx = np.r_[tms[0]/np.sqrt(tms[1]/tms[0]), np.sqrt(tms[:-1]*tms[1:]),
                   tms[-1]*np.sqrt(tms[-1]/tms[-2])]
        my = np.r_[mcs[0]-(mcs[1]-mcs[0])/2, (mcs[:-1]+mcs[1:])/2,
                   mcs[-1]+(mcs[-1]-mcs[-2])/2]
        shape = (len(mcs), len(tms))
        values = np.full((3, *shape), np.nan)
        overlay = np.full(shape, np.nan)
        missing = []
        for g in groups:
            if g['L'] != L: continue
            j = list(tms).index(g['tm_over_tc']); i = list(mcs).index(g['mc'])
            for arm, key in enumerate(('chi0_means', 'chi1_means')):
                if len(g[key]) == 2: values[arm, i, j] = np.mean(g[key])
            if g['unique_chi'] is not None: values[2, i, j] = g['unique_chi']
            if g['classification'] == 'initialization_dependent': overlay[i, j] = 0
            elif g['classification'] in ('unresolved', 'seed_sensitive'): overlay[i, j] = 1
            elif g['classification'] == 'missing': missing.append((tms[j], mcs[i]))
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), sharex=True, sharey=True)
        titles = [r'$\chi_0=0$: late mean', r'$\chi_0=1$: late mean', 'Joint finite-time classification']
        for k, (ax, title) in enumerate(zip(axes, titles)):
            ax.pcolormesh(tx, my, values[k], cmap='bwr', vmin=0, vmax=1,
                          shading='flat', rasterized=True)
            if k == 2:
                ax.pcolormesh(tx, my, overlay, cmap=ListedColormap(['#8160ab', '#bdbdbd']),
                              vmin=0, vmax=1, shading='flat', rasterized=True)
            xx, yy = np.meshgrid(tms, mcs)
            ax.scatter(xx.ravel(), yy.ravel(), s=3, color='.25', alpha=.4)
            if missing:
                a = np.array(missing); ax.scatter(a[:, 0], a[:, 1], marker='x', s=14, color='.2')
            ax.set(xscale='log', xlim=(tms[0], tms[-1]), ylim=(mcs[0], mcs[-1]),
                   xlabel=r'$\tau_m/\tau_c$', title=title)
            ticks = [v for v in [0.3, 1, 3, 10, 30, 50] if min(tms) <= v <= max(tms)]
            ax.set_xticks(ticks, labels=[f'{v:g}' for v in ticks])
        axes[0].set_ylabel(r'$m_c$')
        axes[2].legend(handles=[Patch(facecolor='#8160ab', label='Initial-state dependence'),
                                Patch(facecolor='#bdbdbd', label='Unresolved / seed-sensitive')],
                       loc='best', fontsize=7, framealpha=.9)
        fig.suptitle(f'Empirical phase map: L={L}, two seeds per initialization', fontsize=13)
        fig.subplots_adjust(left=.065, right=.98, top=.83, bottom=.27, wspace=.12)
        ca = fig.add_axes([.29, .14, .42, .025])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap='bwr'),
                          cax=ca, orientation='horizontal')
        cb.set_ticks([0, .5, 1], labels=['0: high activity', '0.5: intermediate', '1: low activity'])
        fig.text(.065, .035, 'After release: 2000–3000 reference tau_c. Tail: 500–1000 tau_c. '
                 'pmem=0.016838; tau_chi=0.3 tau_c.\n'
                 'Dots: sampled points; cells are a display convention. x: missing cases. '
                 'Mean panels retain flagged cases; classification uses convergence diagnostics.', fontsize=8)
        for ext in ('png', 'pdf'): fig.savefig(out/f'phase_{Lname}.{ext}', dpi=200)
        plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('out', type=Path)
    ap.add_argument('--summary', action='store_true')
    ap.add_argument('--manifest', type=Path)
    args = ap.parse_args()
    if args.summary:
        if args.manifest is None: ap.error('--summary needs --manifest')
        aggregate(args.input, args.out, args.manifest)
    else:
        reduce_case(args.input, args.out)
