#!/usr/bin/env python3
"""tau_chi scan at mc=0.2287, tau_m=0.3 tau_c: per-run reduction, summary and movies.

  cw_tchi_scan.py <case_in> <case_out>                       per-run reduction
  cw_tchi_scan.py <results_root> <out> --summary --manifest M campaign summary
  cw_tchi_scan.py <case_in> <movie_out> --render --summary S  native u/P/m/chi movie

The per-run reduction is the four-start reducer (terminal <chi>, drift flags, release
state) plus the phenotype clock read from the case name and checked against the run.
The movie is the four-start native renderer, tagged with tau_chi.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import cw_four_init_scan as four
from cw_mc_tm_scan import TC, parameters

COLORS = {'0': '#63B8C6', '1': '#E8A0A6'}
LABELS = {'0': r'uniform $\chi=0$', '1': r'uniform $\chi=1$'}


def tchi_of(name):
    match = re.search(r'_tc([0-9]+(?:p[0-9]+)?)_tm', name)
    if match is None:
        raise ValueError('Missing tau_chi coordinate in case name')
    return float(match.group(1).replace('p', '.'))


def reduce_case(root, out):
    four.reduce_case(root, out)
    tc = tchi_of(root.name)
    p = parameters(root)
    if not np.isclose(float(p['tau_chi'])/TC, tc, rtol=5e-6):
        raise ValueError('Case path and serialized tau_chi disagree')
    row = json.loads((out/'phase_case.json').read_text())
    s = np.load(out/'series.npz')
    t, chi, std = s['observation_time_tc'], s['chi_mean'], s['chi_std']
    tail = (t >= row['observation_tc']-row['window_tc']) & (t <= row['observation_tc'])
    var = std[tail]**2
    bern = chi[tail]*(1-chi[tail])
    row.update(tau_chi_over_tc=tc, tau_chi_serialized_over_tc=float(p['tau_chi'])/TC,
               tail_chi_spatial_std=float(np.mean(std[tail])),
               # Var_x(chi)/chibar(1-chibar): 1 for a binary field, 0 for a uniform one.
               tail_binariness=float(np.mean(var)/max(np.mean(bern), 1e-12)),
               tchi_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    s.close()
    (out/'phase_case.json').write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')


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
            if not np.isclose(row['tau_chi_over_tc'], item['tau_chi_over_tc'], rtol=1e-10):
                raise ValueError('Wrong tau_chi')
            if not np.isclose(row['tm_over_tc'], item['tm_over_tc'], rtol=1e-10):
                raise ValueError('Wrong memory time')
            rows.append(row)
        except (KeyError, ValueError) as exc:
            invalid.append({'case': item['case'], 'error': str(exc)})
    groups = []
    for tc in manifest['tchi_grid']:
        rr = {r['initialization']: r for r in rows if np.isclose(r['tau_chi_over_tc'], tc)}
        g = {'tau_chi_over_tc': tc,
             'tail_means': {s: rr[s]['tail_mean'] for s in rr},
             'binariness': {s: rr[s]['tail_binariness'] for s in rr},
             'flags': {s: rr[s]['diagnostic_flags'] for s in rr if rr[s]['diagnostic_flags']}}
        if len(rr) == 2:
            g['start_gap'] = rr['1']['tail_mean']-rr['0']['tail_mean']
        groups.append(g)
    summary = {'campaign': manifest['campaign'], 'expected_cases': len(manifest['cases']),
               'available_cases': len(rows), 'missing': missing, 'invalid': invalid,
               'groups': groups, 'cases': rows, 'protocol': manifest['protocol']}
    out.mkdir(parents=True, exist_ok=True)
    (out/'tchi_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    if rows:
        make_figure(rows, manifest, root, out)
    (out/'README.md').write_text(
        '# tau_chi scan at mc=0.2287, tau_m=0.3 tau_c\n\n'
        f'Available {len(rows)}/{len(manifest["cases"])}; missing {len(missing)}, '
        f'invalid {len(invalid)}. One seed per start.\n\n'
        'tchi_scan.png: (a) mean of the exact full-domain <chi> over the last 500 of 1000 '
        'tau_c after release; (b) binariness Var_x(chi)/[chibar(1-chibar)] over the same window '
        '(1 = two-valued field, 0 = uniform); (c) <chi>(t) for every run. Open symbols carry '
        'diagnostic flags (drift / block variation), listed in tchi_summary.json.\n')
    print(json.dumps({'available': len(rows), 'expected': len(manifest['cases']),
                      'missing': len(missing), 'invalid': len(invalid)}))


def make_figure(rows, manifest, root, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    grid = np.array(manifest['tchi_grid'])
    for start, marker in zip(('0', '1'), ('o', 's')):
        rr = sorted((r for r in rows if r['initialization'] == start),
                    key=lambda r: r['tau_chi_over_tc'])
        x = [r['tau_chi_over_tc'] for r in rr]
        for a, key in zip(ax[:2], ('tail_mean', 'tail_binariness')):
            a.plot(x, [r[key] for r in rr], color=COLORS[start], marker=marker, ms=5,
                   lw=1.5, label=LABELS[start])
            for r in rr:
                if r['diagnostic_flags']:
                    a.scatter(r['tau_chi_over_tc'], r[key], marker=marker, s=60,
                              facecolors='none', edgecolors='#30363B', lw=.8, zorder=5)
    cmap = plt.get_cmap('viridis')
    for r in rows:
        s = np.load(root/r['case']/'series.npz')
        k = int(np.argmin(np.abs(grid-r['tau_chi_over_tc'])))
        ax[2].plot(s['observation_time_tc'], s['chi_mean'], lw=.8,
                   color=cmap(k/max(len(grid)-1, 1)), ls='-' if r['initialization'] == '0' else '--')
        s.close()
    for a in ax[:2]:
        a.set_xscale('log')
        a.set_xlabel(r'$\tau_\chi/\tau_c$')
    ax[0].set(ylabel=r'late-time $\langle\chi\rangle$', ylim=(-.025, 1.025))
    ax[1].set(ylabel=r'Var$_x\chi\,/\,[\bar\chi(1-\bar\chi)]$', ylim=(-.025, 1.025))
    ax[0].legend(frameon=False, fontsize=9)
    ax[2].set(xlabel=r'$t/\tau_c$ after release', ylabel=r'$\langle\chi\rangle$', ylim=(-.025, 1.025))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, len(grid)-1))
    cb = fig.colorbar(sm, ax=ax[2], ticks=range(len(grid)), pad=.02)
    cb.ax.set_yticklabels([f'{g:g}' for g in grid])
    cb.set_label(r'$\tau_\chi/\tau_c$ (solid: $\chi_0=0$, dashed: $\chi_0=1$)')
    fig.suptitle(r'$m_c=0.2287$, $\tau_m=0.3\,\tau_c$, $L=256$, one seed', fontsize=11)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'tchi_scan.{ext}', dpi=220, facecolor='white')
    plt.close(fig)


def render(src, out, summary_path):
    import cw_four_init_board as board
    board.render(src, out, summary_path)
    row = next(r for r in json.loads(summary_path.read_text())['cases']
               if Path(r['case']).name == src.name)
    info = json.loads((out/'video.json').read_text())
    info.update(tau_chi_over_tc=row['tau_chi_over_tc'],
                tail_binariness=row['tail_binariness'],
                tchi_render_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'video.json').write_text(json.dumps(info, allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--summary', nargs='?', const=True, default=None,
                        help='campaign summary (no value), or the summary JSON for --render')
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--render', action='store_true')
    args = parser.parse_args()
    if args.render:
        if not isinstance(args.summary, str):
            parser.error('--render requires --summary <tchi_summary.json>')
        render(args.input, args.out, Path(args.summary))
    elif args.summary:
        if args.manifest is None:
            parser.error('--summary requires --manifest')
        aggregate(args.input, args.out, args.manifest)
    else:
        reduce_case(args.input, args.out)
