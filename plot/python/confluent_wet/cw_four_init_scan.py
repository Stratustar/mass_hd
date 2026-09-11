#!/usr/bin/env python3
"""Long-time outcomes of all0/all1/stripe/noise preparations at fixed mc.

Reuse the phase-scan terminal diagnostics, while retaining mixed-start labels,
their float chi0, and their measured state at feedback release.
"""
import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import re

import numpy as np
import cw_mc_tm_scan as base

STARTS = ('0', '1', 'stripe01', 'noise')
LABELS = {'0': 'All 0', '1': 'All 1', 'stripe01': 'Left 0 / right 1', 'noise': 'Gaussian noise'}
COLORS = {'0': '#63B8C6', '1': '#E8A0A6', 'stripe01': '#8C6BB1', 'noise': '#596879'}


def reduce_case(root, out):
    match = re.search(r'_chi(0|1|stripe01|noise)_rep[12]$', root.name)
    if not match:
        raise ValueError('Missing four-start initialization label')
    start = match.group(1)
    p = base.parameters(root)
    config = 'stripe' if start == 'stripe01' else 'noise' if start == 'noise' else 'uniform'
    if p['chi_config'] != config:
        raise ValueError('Initialization label disagrees with serialized chi_config')
    expected_chi0 = float(start) if start in ('0', '1') else .5
    if float(p['chi0']) != expected_chi0:
        raise ValueError('Initialization chi0 disagrees with case label')
    # The legacy reducer identifies tm from the shared _tm..._chi naming scheme.
    # It casts chi0 to int, so explicitly restore mixed starts before publishing.
    with redirect_stdout(io.StringIO()):
        base.reduce_case(root, out)
    row = json.loads((out/'phase_case.json').read_text())
    series = np.load(out/'series.npz')
    freeze = row['preparation_steps']
    before = np.flatnonzero(series['t'] <= freeze)
    if not len(before):
        raise ValueError('No pre-release sample')
    j = int(before[-1])
    row.update(initialization=start, chi0=expected_chi0,
        chi_config=config, chi_seed=int(p['chi_seed']),
        base_analysis_sha256=row['script_sha256'],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        initial_sample_step=float(series['t'][0]),
        initial_chi_mean=float(series['chi_mean'][0]),
        initial_chi_std=float(series['chi_std'][0]),
        release_sample_step=float(series['t'][j]),
        release_sample_offset_tc=float((series['t'][j]-freeze)/base.TC),
        release_chi_mean=float(series['chi_mean'][j]),
        release_chi_std=float(series['chi_std'][j]),
        release_chi_variance=float(series['chi_std'][j]**2),
        release_m_mean=float(series['m_mean'][j]),
        preparation_note='Only chi reaction is frozen. Mixed t=0 patterns advect/diffuse before release.')
    series.close()
    (out/'phase_case.json').write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'case': row['case'], 'initialization': start,
        'release_chi_std': row['release_chi_std'], 'tail_mean': row['tail_mean'],
        'flags': row['diagnostic_flags']}))


def classify_four(rows):
    if len(rows) != 8:
        return 'missing'
    if any(not r['settled_by_diagnostics'] for r in rows):
        return 'unresolved'
    arms = [np.array([r['tail_mean'] for r in sorted(rows, key=lambda r: r['replicate'])
                      if r['initialization'] == start]) for start in STARTS]
    if any(len(a) != 2 for a in arms):
        return 'missing'
    if any(np.ptp(a) > .08 for a in arms):
        return 'seed_sensitive'
    if np.ptp(np.concatenate(arms)) <= .05:
        return 'initialization_independent'
    for i, a in enumerate(arms):
        for b in arms[i+1:]:
            lo, hi = (a, b) if np.mean(a) < np.mean(b) else (b, a)
            if np.min(hi-lo) > .20 and np.min(hi)-np.max(lo) > .10:
                return 'initialization_dependent'
    return 'unresolved'


def aggregate(root, out, manifest_path):
    manifest = json.loads(manifest_path.read_text())
    rows, missing, invalid = [], [], []
    for item in manifest['cases']:
        path = root/item['case']/'phase_case.json'
        if not path.exists():
            missing.append(item['case'])
            continue
        try:
            row = json.loads(path.read_text())
            for key in ('case', 'L', 'mc', 'initialization', 'chi0', 'replicate', 'seed',
                        'chi_seed', 'nsteps', 'preparation_steps'):
                if row[key] != item[key]:
                    raise ValueError(f'Wrong {key}')
            if not np.isclose(row['tm_over_tc'], item['tm_over_tc'], rtol=1e-10):
                raise ValueError('Wrong memory time')
            for key in ('pmem', 'tau_chi'):
                if row[key] != manifest['fixed'][key]:
                    raise ValueError(f'Wrong {key}')
            rows.append(row)
        except (KeyError, ValueError) as exc:
            invalid.append({'case': item['case'], 'error': str(exc)})
    groups = []
    for tm in manifest['tm_grid']:
        rr = [r for r in rows if np.isclose(r['tm_over_tc'], tm, rtol=1e-10)]
        arms = {s: [r for r in rr if r['initialization'] == s] for s in STARTS}
        uniform = arms['0']+arms['1']
        groups.append({'tm_over_tc': tm, 'mc': manifest['fixed']['mc'],
            'four_start_classification': classify_four(rr),
            'uniform_pair_classification': base.classify(uniform)[0],
            'arms': {s: {'tail_means': [r['tail_mean'] for r in arms[s]],
                        'release_chi_means': [r['release_chi_mean'] for r in arms[s]],
                        'release_chi_stds': [r['release_chi_std'] for r in arms[s]],
                        'cases': [r['case'] for r in arms[s]]} for s in STARTS},
            'flags': {r['case']: r['diagnostic_flags'] for r in rr if r['diagnostic_flags']}})
    summary = {'campaign': manifest['campaign'], 'expected_cases': len(manifest['cases']),
        'available_cases': len(rows), 'missing': missing, 'invalid': invalid, 'groups': groups,
        'cases': rows, 'protocol': manifest['protocol'], 'criteria': manifest['statistics'],
        'interpretation': manifest['interpretation']}
    out.mkdir(parents=True, exist_ok=True)
    (out/'four_init_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    make_figure(rows, manifest, out)
    (out/'README.md').write_text(
        '# Four initializations at mc=0.2287\n\n'
        f'Available {len(rows)}/{len(manifest["cases"])}; missing {len(missing)}, invalid {len(invalid)}.\n\n'
        'Curves show each initialization separately: mean of two seed outcomes, with shading '
        'for their min-to-max range (not confidence intervals). Open symbols mark flagged cases. '
        'All outcomes use the last 500 reference tau_c of the 2000-tau_c observation period.\n\n'
        'Preparation follows the phase scan: chi reaction off, transport on. Left/right and '
        'noise describe t=0; use release_chi_mean and release_chi_std to inspect how much '
        'mixing has already occurred before feedback release.\n\n'
        'JSON retains flags, individual seeds and both the uniform-pair and four-start '
        'finite-time classifications. Missing or unsettled cases are never silently discarded.\n')
    print(json.dumps({'available': len(rows), 'expected': len(manifest['cases']),
                      'missing': len(missing), 'invalid': len(invalid)}))


def make_figure(rows, manifest, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
        'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    tm = np.array(manifest['tm_grid'])
    for start, marker in zip(STARTS, ('o', 's', '^', 'D')):
        lo, hi, mean = [], [], []
        for t in tm:
            arm = [r['tail_mean'] for r in rows if r['initialization'] == start
                   and np.isclose(r['tm_over_tc'], t, rtol=1e-10)]
            lo.append(min(arm) if arm else np.nan)
            hi.append(max(arm) if arm else np.nan)
            mean.append(float(np.mean(arm)) if len(arm) == 2 else np.nan)
        ax.fill_between(tm, lo, hi, color=COLORS[start], alpha=.12, linewidth=0)
        ax.plot(tm, mean, color=COLORS[start], marker=marker, ms=4, lw=1.5, label=LABELS[start])
        for r in rows:
            if r['initialization'] == start and not r['settled_by_diagnostics']:
                ax.scatter(r['tm_over_tc'], r['tail_mean'], marker=marker, s=48,
                           facecolors='none', edgecolors='#30363B', linewidths=.8, zorder=5)
    ax.set(xlabel=r'$\tau_m/\tau_c$', ylabel=r'Late-time $\langle\chi\rangle$',
           xlim=(-.05, 18.35), ylim=(-.025, 1.025), title=r'$m_c=0.2287$, $L=256$, two seeds per start')
    ax.legend(frameon=False, ncol=2, fontsize=9)
    fig.subplots_adjust(bottom=.22, left=.12, right=.97, top=.88)
    fig.text(.12, .06, 'Last 500 of 2000 reference tau_c after preparation. Shading: two-seed range.\n'
             'Open symbols: diagnostic flags. Mixed starts are advected during preparation.', fontsize=8, color='#596879')
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'four_initializations.{ext}', dpi=220, facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--summary', action='store_true')
    parser.add_argument('--manifest', type=Path)
    args = parser.parse_args()
    if args.summary:
        if args.manifest is None:
            parser.error('--summary requires --manifest')
        aggregate(args.input, args.out, args.manifest)
    else:
        reduce_case(args.input, args.out)
