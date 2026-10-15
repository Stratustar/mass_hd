#!/usr/bin/env python3
"""Seed-averaged summary of the hill-law family scans (20261015 nscan and tchiscan).

  cw_hill_family_scan.py <out> --scan <results_root> <manifest> [--scan ...]
                         [--legacy four_init_summary.json] [--cell cells_hill_*.json]

Per-run reduction is cw_hill_tm_scan.py (phase_case.json). Here each (member, tau_m, start)
is averaged over its seeds: mean of the terminal <chi>, seed min/max, binariness, Var_t<chi>.
Figures (hill convention, chi = 1 active; smooth = solid, step = dashed, single cell = points):
  nscan_chi    one curve per Hill exponent n; the step limit (h-n = 0, wp = 0) dashed black and
               the legacy step-law scan (rep 1-2, 1 - <chi>) dashed grey
  tchiscan_chi one curve per tau_chi (0.3 tau_c taken from the nscan n = 16 member) and the
               single-cell model as points
Missing runs are listed, never silently dropped; a group with fewer seeds says so.
"""
import argparse
import json
from pathlib import Path

import numpy as np

STARTS = ('1', '0')


def load(scans):
    groups, missing, members = {}, [], {}
    for root, manifest_path in scans:
        manifest = json.loads(Path(manifest_path).read_text())
        for item in manifest['cases']:
            members[item['member']] = {k: item[k] for k in ('h_n', 'wp', 'tau_chi_over_tc')}
            path = Path(root)/item['case']/'phase_case.json'
            if not path.exists():
                missing.append(item['case'])
                continue
            r = json.loads(path.read_text())
            if r['case'] != item['case'] or r['replicate'] != item['replicate'] or \
                    not np.isclose(r['tm_over_tc'], item['tm_over_tc'], rtol=1e-9):
                raise ValueError(f'Row does not match manifest: {path}')
            key = (item['member'], item['tm_over_tc'], item['initialization'])
            groups.setdefault(key, []).append(r)
    return groups, missing, members


def summarise(groups, members):
    out = []
    for (member, tm, start), rr in sorted(groups.items()):
        x = np.array([r['tail_mean'] for r in rr])
        out.append({'member': member, **members[member], 'tm_over_tc': tm, 'initialization': start,
                    'n_seeds': len(rr), 'chi_mean': float(x.mean()), 'chi_min': float(x.min()),
                    'chi_max': float(x.max()), 'chi_seeds': x.tolist(),
                    'binariness_mean': float(np.mean([r['tail_binariness'] for r in rr])),
                    'var_t_mean': float(np.mean([r['tail_std_time']**2 for r in rr])),
                    'flags': {r['case']: r['diagnostic_flags'] for r in rr if r['diagnostic_flags']}})
    return out


def curve(rows, member, start):
    rr = sorted((r for r in rows if r['member'] == member and r['initialization'] == start),
                key=lambda r: r['tm_over_tc'])
    return (np.array([r['tm_over_tc'] for r in rr]), np.array([r['chi_mean'] for r in rr]),
            np.array([r['chi_min'] for r in rr]), np.array([r['chi_max'] for r in rr]))


def legacy_curve(path):
    """Legacy step scan, rep 1-2 averaged, as active fraction per hill start label."""
    if not path or not Path(path).exists():
        return {}
    acc = {}
    for r in json.loads(Path(path).read_text())['cases']:
        if r['initialization'] not in ('0', '1'):
            continue
        start = '1' if r['initialization'] == '0' else '0'
        acc.setdefault((start, r['tm_over_tc']), []).append(1 - r['tail_mean'])
    out = {}
    for (start, tm), v in acc.items():
        out.setdefault(start, []).append((tm, float(np.mean(v))))
    return {s: np.array(sorted(v)) for s, v in out.items()}


def cell_curve(path):
    if not path or not Path(path).exists():
        return {}
    out = {}
    for r in json.loads(Path(path).read_text())['runs']:
        out.setdefault('1' if r['start'] == 'active' else '0', []).append((r['tau_m'], r['tail_mean']))
    return {s: np.array(sorted(v)) for s, v in out.items()}


def figures(rows, out, legacy, cell):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.labelsize': 12,
                         'legend.fontsize': 8, 'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
    present = {r['member'] for r in rows}

    def draw(ax, member, color, label):
        for start, mk, mfc in (('1', 'o', color), ('0', 's', 'white')):
            t, m, lo, hi = curve(rows, member, start)
            if len(t):
                ax.errorbar(t, m, yerr=[m - lo, hi - m], fmt='-' + mk, color=color, mfc=mfc,
                            mec=color, ms=4, lw=1.4, elinewidth=.8, capsize=0,
                            label=label if start == '1' else None)

    def axes(ax):
        ax.set_xlabel(r'$\tau_m/\tau_c$'); ax.set_ylabel(r'$\langle\chi\rangle$')
        ax.set_xlim(-.3, 18.8); ax.set_ylim(-.03, 1.03)

    # --- n scan
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ns = sorted({int(m.split('_')[0][2:]) for m in present if m.startswith('hn') and m.endswith('tc0p3')
                 and not m.startswith('hn0_')})
    cmap = plt.get_cmap('viridis')
    for k, n in enumerate(ns):
        draw(ax, f'hn{n}_wp0p002_tc0p3', cmap(.1 + .8*k/max(len(ns) - 1, 1)), f'n = {n}')
    for start in ('1', '0'):
        t, m, _, _ = curve(rows, 'hn0_wp0_tc0p3', start)
        if len(t):
            ax.plot(t, m, '--', color='#30363B', lw=1.5, label='step limit' if start == '1' else None)
        if start in legacy:
            ax.plot(legacy[start][:, 0], legacy[start][:, 1], '--', color='#9AA3AB', lw=1.2,
                    label='legacy step' if start == '1' else None)
    axes(ax); ax.legend(loc='center right', bbox_to_anchor=(1, .5))
    fig.tight_layout(pad=.3)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out/f'nscan_chi.{ext}', dpi=300)
    plt.close(fig)

    # --- tau_chi scan (0.3 tau_c = the nscan n = 16 member)
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    tcs = sorted({m.split('_tc')[1] for m in present if m.startswith('hn16_wp0p002_tc')},
                 key=lambda s: float(s.replace('p', '.')))
    cmap = plt.get_cmap('plasma')
    for k, tc in enumerate(tcs):
        draw(ax, f'hn16_wp0p002_tc{tc}', cmap(.05 + .8*k/max(len(tcs) - 1, 1)),
             rf'$\tau_\chi={tc.replace("p", ".")}\,\tau_c$')
    for start, mfc in (('1', '#21918c'), ('0', 'white')):
        if start in cell:
            ax.plot(cell[start][:, 0], cell[start][:, 1], ls='none', marker='D', ms=4, color='#21918c',
                    mfc=mfc, mew=1.1, label='single cell' if start == '1' else None, zorder=5)
    axes(ax); ax.legend(loc='center right', bbox_to_anchor=(1, .5))
    fig.tight_layout(pad=.3)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out/f'tchiscan_chi.{ext}', dpi=300)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('out', type=Path)
    p.add_argument('--scan', nargs=2, action='append', required=True, metavar=('RESULTS_ROOT', 'MANIFEST'))
    p.add_argument('--legacy', type=Path)
    p.add_argument('--cell', type=Path)
    a = p.parse_args()
    groups, missing, members = load(a.scan)
    rows = summarise(groups, members)
    a.out.mkdir(parents=True, exist_ok=True)
    short = [r for r in rows if r['n_seeds'] < 3]
    (a.out/'family_summary.json').write_text(json.dumps(
        {'scans': [[str(x), str(y)] for x, y in a.scan], 'missing': missing,
         'groups_with_fewer_seeds': [(r['member'], r['tm_over_tc'], r['initialization'], r['n_seeds']) for r in short],
         'legacy_reference': str(a.legacy) if a.legacy else None,
         'cell_reference': str(a.cell) if a.cell else None, 'rows': rows}, indent=1))
    if rows:
        figures(rows, a.out, legacy_curve(a.legacy), cell_curve(a.cell))
    print(json.dumps({'groups': len(rows), 'missing_runs': len(missing), 'short_groups': len(short)}))


if __name__ == '__main__':
    main()
