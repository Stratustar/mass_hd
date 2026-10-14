#!/usr/bin/env python3
"""Hill response law, tau_m scan at mc=0.2287: per-run reduction, summary and movies.

  cw_hill_tm_scan.py <case_in> <case_out>                        per-run reduction
  cw_hill_tm_scan.py <results_root> <out> --summary --manifest M [--reference R]
                                                                 campaign summary
  cw_hill_tm_scan.py <case_in> <movie_out> --render --summary S  native u/P/m/chi movie
  cw_hill_tm_scan.py <fetched_root> <page_dir> --board --summary S [--reference R]
                                                                 offline video board

The reduction is the four-start reducer (terminal <chi>, drift flags, release state). The
hill law FLIPS the phenotype convention -- chi = 1 is active -- so <chi> here is the ACTIVE
fraction, the legacy reducer's chi < 0.2 / chi > 0.8 occupancies are swapped back into
their physical names, and the optional legacy reference R (the 20260911 four_init_summary)
is drawn as 1 - <chi_legacy> against the same tau_m with its starts relabelled.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import cw_four_init_scan as four
from cw_mc_tm_scan import TC, parameters

COLORS = {'1': '#E07A5F', '0': '#3D7EA6'}
LABELS = {'1': r'uniform $\chi=1$ (active)', '0': r'uniform $\chi=0$ (passive)'}
LEGACY = {'0': '1', '1': '0'}      # legacy start label -> hill start label


def reduce_case(root, out):
    p = parameters(root)
    if p.get('switch_model') != 'hill':
        raise ValueError('Not a switch-model=hill run')
    four.reduce_case(root, out)
    row = json.loads((out/'phase_case.json').read_text())
    s = np.load(out/'series.npz')
    t, chi, std, m = s['observation_time_tc'], s['chi_mean'], s['chi_std'], s['m_mean']
    tail = (t >= row['observation_tc']-row['window_tc']) & (t <= row['observation_tc'])
    var, bern = std[tail]**2, chi[tail]*(1-chi[tail])
    # The legacy reducer labels chi < 0.2 'active'; under the hill law that is passive.
    active, passive = row['tail_passive_occupancy'], row['tail_active_occupancy']
    row.update(chi_convention='hill: chi=1 active', tail_active_fraction=row['tail_mean'],
               tail_active_occupancy=active, tail_passive_occupancy=passive,
               tail_m=float(np.mean(m[tail])),
               tail_binariness=float(np.mean(var)/max(np.mean(bern), 1e-12)),
               hill={k: p[k] for k in ('h_n', 'wp', 'beta', 'g_n', 'g_chic', 'g_lock',
                                       'Dm', 'Dchi', 'Dbio')},
               hill_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    s.close()
    (out/'phase_case.json').write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')


def tchi_text(fixed, tex=False):
    """The phenotype clock of a campaign: a fixed number of steps, or a multiple of tau_m."""
    if 'tau_chi_over_tm' in fixed:
        r = fixed['tau_chi_over_tm']
        return rf'$\tau_\chi={r:g}\,\tau_m$' if tex else f'{r:g} τm'
    x = fixed.get('tau_chi_steps', 202.3)/TC
    return rf'$\tau_\chi={x:.2g}\,\tau_c$' if tex else f'{x:.2g} τc'


def classify(rr):
    """One seed per start: dependent if the active start ends > 0.2 above the passive one."""
    if set(rr) != {'0', '1'}:
        return 'missing'
    if not all(r['settled_by_diagnostics'] for r in rr.values()):
        return 'unresolved'
    gap = rr['1']['tail_mean'] - rr['0']['tail_mean']
    return ('initialization_dependent' if gap > .2 else
            'initialization_independent' if abs(gap) <= .05 else 'unresolved')


def legacy_rows(path):
    if path is None:
        return []
    rows = [r for r in json.loads(Path(path).read_text())['cases']
            if r['replicate'] == 1 and r['initialization'] in LEGACY]
    return [{'tm_over_tc': r['tm_over_tc'], 'initialization': LEGACY[r['initialization']],
             'tail_active_fraction': 1-r['tail_mean'], 'tail_m': r['tail_m_mean'],
             'case': r['case']} for r in rows]


def aggregate(root, out, manifest_path, reference=None, compare=None):
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
            for key in ('h-n', 'wp', 'beta'):
                if row['hill'][key.replace('-', '_')] != manifest['fixed'][key]:
                    raise ValueError(f'Wrong {key}')
            rows.append(row)
        except (KeyError, ValueError) as exc:
            invalid.append({'case': item['case'], 'error': str(exc)})
    legacy = legacy_rows(reference)
    # Another hill campaign drawn for comparison (e.g. a different tau_chi), as its own rows.
    other = json.loads(Path(compare).read_text()) if compare else None
    groups = []
    for tm in manifest['tm_grid']:
        rr = {r['initialization']: r for r in rows if np.isclose(r['tm_over_tc'], tm, rtol=1e-10)}
        lg = {r['initialization']: r for r in legacy if np.isclose(r['tm_over_tc'], tm, rtol=1e-9)}
        g = {'tm_over_tc': tm, 'classification': classify(rr),
             'tail_active_fraction': {s: rr[s]['tail_mean'] for s in rr},
             'tail_m': {s: rr[s]['tail_m'] for s in rr},
             'binariness': {s: rr[s]['tail_binariness'] for s in rr},
             'legacy_active_fraction': {s: lg[s]['tail_active_fraction'] for s in lg},
             'flags': {s: rr[s]['diagnostic_flags'] for s in rr if rr[s]['diagnostic_flags']}}
        if len(rr) == 2:
            g['start_gap'] = rr['1']['tail_mean'] - rr['0']['tail_mean']
        groups.append(g)
    summary = {'campaign': manifest['campaign'], 'expected_cases': len(manifest['cases']),
               'available_cases': len(rows), 'missing': missing, 'invalid': invalid,
               'convention': manifest['convention'], 'fixed': manifest['fixed'],
               'legacy_reference': str(reference) if reference else None,
               'compare': str(compare) if compare else None,
               'groups': groups, 'cases': rows, 'protocol': manifest['protocol']}
    out.mkdir(parents=True, exist_ok=True)
    (out/'hill_tm_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    if rows:
        make_figure(rows, legacy, manifest, root, out, other)
    (out/'README.md').write_text(
        '# Hill law (beta = 0) tau_m scan at mc = 0.2287\n\n'
        f'Available {len(rows)}/{len(manifest["cases"])}; missing {len(missing)}, '
        f'invalid {len(invalid)}. One seed per start. chi = 1 is ACTIVE.\n\n'
        'hill_tm_scan.png: (a) late-time active fraction <chi> (window max(500, 20 tau_m) '
        'tau_c at the end of the 2000 tau_c observation) against tau_m, both starts; grey: '
        'the legacy step-law partner runs, 1 - <chi_legacy>; (b) late <m>, same layout; '
        '(c) <chi>(t) of every run. Open symbols carry diagnostic flags.\n')
    print(json.dumps({'available': len(rows), 'expected': len(manifest['cases']),
                      'missing': len(missing), 'invalid': len(invalid)}))


def make_figure(rows, legacy, manifest, root, out, other=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for start, marker in (('1', 'o'), ('0', 's')):
        for a, key, lkey in ((ax[0], 'tail_mean', 'tail_active_fraction'),
                             (ax[1], 'tail_m', 'tail_m')):
            lg = sorted((r for r in legacy if r['initialization'] == start),
                        key=lambda r: r['tm_over_tc'])
            if lg:
                a.plot([r['tm_over_tc'] for r in lg], [r[lkey] for r in lg], color='#9AA3AB',
                       marker=marker, ms=3.5, lw=1, ls='--', zorder=1,
                       label='legacy step law' if start == '1' else None)
            if other:
                oo = sorted((r for r in other['cases'] if r['initialization'] == start),
                            key=lambda r: r['tm_over_tc'])
                a.plot([r['tm_over_tc'] for r in oo], [r[key] for r in oo], color=COLORS[start],
                       marker=marker, ms=3.5, lw=1, ls=':', alpha=.55, zorder=2,
                       label=f'hill, τχ = {tchi_text(other["fixed"])}' if start == '1' else None)
            rr = sorted((r for r in rows if r['initialization'] == start),
                        key=lambda r: r['tm_over_tc'])
            a.plot([r['tm_over_tc'] for r in rr], [r[key] for r in rr], color=COLORS[start],
                   marker=marker, ms=5, lw=1.6, label=LABELS[start], zorder=3)
            for r in rr:
                if r['diagnostic_flags']:
                    a.scatter(r['tm_over_tc'], r[key], marker=marker, s=70, facecolors='none',
                              edgecolors='#30363B', lw=.8, zorder=5)
    ax[1].axhline(manifest['fixed']['mc'], color='#30363B', lw=.7, ls=':')
    grid = np.array(manifest['tm_grid'])
    cmap = plt.get_cmap('viridis')
    for r in rows:
        s = np.load(root/r['case']/'series.npz')
        k = int(np.argmin(np.abs(grid-r['tm_over_tc'])))
        ax[2].plot(s['observation_time_tc'], s['chi_mean'], lw=.8,
                   color=cmap(k/max(len(grid)-1, 1)), ls='-' if r['initialization'] == '1' else '--')
        s.close()
    for a in ax[:2]:
        a.set_xlabel(r'$\tau_m/\tau_c$')
    ax[0].set(ylabel=r'late-time active fraction $\langle\chi\rangle$', ylim=(-.025, 1.025))
    ax[1].set(ylabel=r'late-time $\langle m\rangle$  (dotted: $m_c$)')
    ax[0].legend(frameon=False, fontsize=8.5)
    ax[2].set(xlabel=r'$t/\tau_c$ after release', ylabel=r'$\langle\chi\rangle$', ylim=(-.025, 1.025))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, len(grid)-1))
    cb = fig.colorbar(sm, ax=ax[2], ticks=range(len(grid)), pad=.02)
    cb.ax.set_yticklabels([f'{g:.3g}' for g in grid])
    cb.set_label(r'$\tau_m/\tau_c$ (solid: $\chi_0=1$, dashed: $\chi_0=0$)')
    f = manifest['fixed']
    fig.suptitle(rf'hill law: $n_1={f["h-n"]:g}$, $w_p={f["wp"]:g}$, $\beta={f["beta"]:g}$, '
                 + tchi_text(f, tex=True) + r'; $m_c=0.2287$, $L=256$, one seed', fontsize=11)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'hill_tm_scan.{ext}', dpi=220, facecolor='white')
    plt.close(fig)


def render(src, out, summary_path):
    import cw_four_init_board as board
    board.render(src, out, summary_path)
    row = next(r for r in json.loads(summary_path.read_text())['cases']
               if Path(r['case']).name == src.name)
    info = json.loads((out/'video.json').read_text())
    info.update(chi_convention=row['chi_convention'], tail_binariness=row['tail_binariness'],
                hill=row['hill'],
                hill_render_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'video.json').write_text(json.dumps(info, allow_nan=False))


def board(data_root, out, summary_path, reference=None):
    """Offline page over fetched <data_root>/<case>/dashboard/{fields.mp4,poster.png,video.json}.

    The memory times have different preparation lengths (max(100 tau_c, 10 tau_m)), so the
    movies are synchronised on time after release, not on frame index; the page maps that
    clock into each movie separately. reference: the legacy 20260911 four_init_summary,
    drawn on the summary panels as 1 - <chi_legacy> with its starts relabelled.
    """
    import os
    summary = json.loads(Path(summary_path).read_text())
    if summary['missing'] or summary['invalid']:
        raise ValueError('Board needs a complete summary')
    runs = []
    for row in summary['cases']:
        meta = Path(data_root)/row['case']/'dashboard/video.json'
        r = json.loads(meta.read_text())
        for key in ('case', 'initialization', 'tm_over_tc', 'preparation_steps', 'tail_mean'):
            if r[key] != row[key]:
                raise ValueError(f'Wrong {key}: {meta}')
        movie, poster = meta.with_name('fields.mp4'), meta.with_name('poster.png')
        if movie.stat().st_size != r['video_bytes'] or \
                hashlib.sha256(movie.read_bytes()).hexdigest() != r['video_sha256']:
            raise ValueError(f'Movie does not match metadata: {movie}')
        if r['frames'] != row['nsteps']//337+1 or len(r['chi']) != r['frames'] or \
                len(r['memory']) != r['frames'] or r['native_shape'] != [256, 256]:
            raise ValueError(f'Frame count or shape mismatch: {meta}')
        runs.append({k: r[k] for k in ('case', 'initialization', 'tm_over_tc', 'tail_mean',
                                       'frames', 'frame_steps', 'fps', 'duration', 'width',
                                       'height', 'preparation_steps', 'chi', 'memory')}
                    | {k: row[k] for k in ('tail_m', 'tail_binariness', 'observation_tc',
                                           'diagnostic_flags')}
                    | {'url': os.path.relpath(movie, out), 'poster': os.path.relpath(poster, out)})
    if len({(r['tm_over_tc'], r['initialization']) for r in runs}) != len(runs) or \
            len({(r['frame_steps'], r['fps'], r['width'], r['height']) for r in runs}) != 1:
        raise ValueError('Duplicate runs or movies that cannot share one clock')
    payload = {'runs': runs, 'legacy': legacy_rows(reference), 'tau_c': TC,
               'mc': summary['fixed']['mc'], 'shape': [256, 256],
               'law': {k: summary['fixed'][k] for k in ('h-n', 'wp', 'beta')}
                      | {'tau_chi': tchi_text(summary['fixed'])},
               'tc_per_second': runs[0]['frame_steps']*runs[0]['fps']/TC,
               # Geometry of the native four-panel composite written by cw_four_init_board.
               'panel': {'size': 256, 'header': 32, 'chi_x': 3*(256+12)}}
    template = Path(__file__).with_name('cw_hill_board.html.in').read_text()
    if template.count('__DATA__') != 1:
        raise ValueError('Invalid template marker')
    text = template.replace('__DATA__', json.dumps(payload).replace('</', '<\\/'))
    out.mkdir(parents=True, exist_ok=True)
    (out/'index.html').write_text(text)
    print(json.dumps({'board': str(out/'index.html'), 'videos': len(runs),
                      'html_MB': len(text.encode())/1e6}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--summary', nargs='?', const=True, default=None,
                        help='campaign summary (no value), or the summary JSON for --render')
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--compare', type=Path,
                        help='another hill_tm_summary.json drawn dotted on the summary figure')
    parser.add_argument('--reference', type=Path,
                        help='legacy 20260911 four_init_summary.json drawn for comparison')
    parser.add_argument('--render', action='store_true')
    parser.add_argument('--board', action='store_true',
                        help='<fetched results root> <page dir> --board --summary <json> '
                             '[--reference <legacy json>]')
    args = parser.parse_args()
    if args.board:
        if not isinstance(args.summary, str):
            parser.error('--board requires --summary <hill_tm_summary.json>')
        board(args.input.resolve(), args.out.resolve(), Path(args.summary), args.reference)
    elif args.render:
        if not isinstance(args.summary, str):
            parser.error('--render requires --summary <hill_tm_summary.json>')
        render(args.input, args.out, Path(args.summary))
    elif args.summary:
        if args.manifest is None:
            parser.error('--summary requires --manifest')
        aggregate(args.input, args.out, args.manifest, args.reference, args.compare)
    else:
        reduce_case(args.input, args.out)
