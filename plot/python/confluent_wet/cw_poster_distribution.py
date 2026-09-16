#!/usr/bin/env python3
"""Native memory/Beta and pressure/Gaussian diagnostics for the poster campaign.

Per case: ROOT is a native archive, OUT is its results directory.
With --summary: ROOT is the campaign results tree, OUT is a summary directory.
Two streaming passes retain moments, histograms and CDFs, never all raw samples.
The model distributions match the measured first and second moments; no fitting,
centering, clipping, KDE or independent-site uncertainty estimates are used.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import betainc, betaincc, ndtr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[3]
DEFAULT_MANIFEST = REPO / 'cases/20260917/cw_poster_distributions/manifest.json'
VERSION = 1
FIELDS = ('m', 'pressure')
ALIASES = {'model': 'model_name', 'bc': 'BC', 'initial-order': 'init_order'}
REQUIRED = ('LX', 'LY', 'nsteps', 'nstart', 'ninfo', 'seed', 'model_name',
            'open_loop', 'zeta_open', 'zeta', 'zeta0_frac', 'tau_m', 'tau_chi',
            'mc', 'pmem', 'pmem_width', 'chi_width', 'Dbio', 'frame_light')
STYLE = {'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
         'pdf.fonttype': 42, 'ps.fonttype': 42, 'legend.frameon': False}
TRACE_COLUMNS = ['step', 't_over_tc'] + [f'{f}_{k}' for f in FIELDS
    for k in ('mean', 'raw_second', 'variance', 'sigma', 'tail', 'minimum', 'maximum')]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def read_native(path):
    payload = path.read_bytes()
    data = json.loads(payload)['data']
    return {key: value['value'] for key, value in data.items()}, hashlib.sha256(payload).hexdigest()


def native_fields(path, sites):
    data, digest = read_native(path)
    values = {}
    for name in FIELDS:
        if name not in data:
            raise ValueError(f'{path}: native {name} field missing')
        x = np.asarray(data[name], dtype=np.float64).ravel()
        if x.size != sites or not np.isfinite(x).all():
            raise ValueError(f'{path}: invalid {name} shape/finite check ({x.size} != {sites})')
        if name == 'm' and (x.min() < 0 or x.max() > 1):
            raise ValueError(f'{path}: m outside [0,1]; data are never clipped')
        values[name] = x
    return values, digest


def check_parameters(par, row):
    missing = sorted(set(REQUIRED) - set(par))
    if missing:
        raise ValueError(f'Missing serialized parameters: {missing}')
    checked, unserialized = [], []
    for key, expected in row['expected_parameters'].items():
        name = ALIASES.get(key, key.replace('-', '_'))
        if name not in par:
            unserialized.append(key)
            continue
        actual = par[name]
        try:
            ok = bool(np.isclose(float(actual), float(expected), rtol=6e-6, atol=1e-12))
        except (TypeError, ValueError):
            ok = str(actual) == str(expected)
        if not ok:
            raise ValueError(f'Runtime mismatch {key}: {actual!r} != {expected!r}')
        checked.append(key)
    if str(par['model_name']) != 'confluent-wet' or int(par['open_loop']) != 1:
        raise ValueError('This analysis requires open-loop confluent-wet')
    return {'checked': checked, 'not_serialized': unserialized,
            'float_rtol': 6e-6, 'float_atol': 1e-12}


def moment_stats(mean, second):
    variance = second - mean * mean
    if variance < -1e-13 * max(abs(second), 1):
        raise ValueError('Negative variance from native moments')
    variance = max(0., variance)
    return {'mean': float(mean), 'raw_second': float(second),
            'variance': float(variance), 'sigma': float(np.sqrt(variance))}


def closure(field, mean, second):
    stats = moment_stats(mean, second)
    variance = stats['variance']
    result = dict(family='Beta' if field == 'm' else 'Gaussian', **stats)
    if variance == 0:
        return dict(result, status='degenerate', reason='Zero native variance: no continuous PDF')
    if field == 'm':
        if not 0 < mean < 1 or variance >= mean * (1 - mean):
            return dict(result, status='undefined', reason='Moments reach/exceed continuous Beta domain')
        concentration = mean * (1 - mean) / variance - 1
        alpha, beta = mean * concentration, (1 - mean) * concentration
        model_mean = alpha / (alpha + beta)
        model_second = alpha * (alpha + 1) / ((alpha + beta) * (alpha + beta + 1))
        result.update(alpha=float(alpha), beta=float(beta), concentration=float(concentration))
    else:
        model_mean, model_second = mean, variance + mean * mean
    residual = [float(model_mean - mean), float(model_second - second)]
    if not np.allclose(residual, 0, atol=1e-12, rtol=0):
        raise ValueError(f'Model does not reproduce native moments: {residual}')
    return dict(result, status='valid', model_mean=float(model_mean),
                model_raw_second=float(model_second), moment_residual=residual)


def model_cdf(model, x):
    x = np.asarray(x, dtype=float)
    if model['status'] != 'valid':
        return np.full(x.shape, np.nan)
    if model['family'] == 'Beta':
        return betainc(model['alpha'], model['beta'], np.clip(x, 0, 1))
    return ndtr((x - model['mean']) / model['sigma'])


def model_bins(model, edges):
    """Bin-integrated probabilities: finite even for alpha/beta below one."""
    cdf = model_cdf(model, edges)
    if model['status'] != 'valid':
        return np.full(len(edges) - 1, np.nan)
    if model['family'] == 'Beta':
        survival = betaincc(model['alpha'], model['beta'], np.clip(edges, 0, 1))
    else:
        survival = ndtr((model['mean'] - edges) / model['sigma'])
    mass = np.where(cdf[1:] < .5, np.diff(cdf), -np.diff(survival))
    if np.any(mass < -1e-14):
        raise ValueError('Negative model bin probability')
    return np.maximum(mass, 0)


def model_tail(model, threshold, field):
    if model['status'] != 'valid':
        return None
    if field == 'm':
        return float(model_cdf(model, threshold))
    return float(ndtr((model['mean'] - threshold) / model['sigma']))


def cdf_counts(x, edges):
    # searchsorted(..., left) counts x == edge in F(edge) = Pr(X <= edge).
    return np.bincount(np.searchsorted(edges, x, side='left'),
                       minlength=len(edges) + 1)[:-1].cumsum()


def distribution_stats(stats, edges, cdf, threshold, field):
    model = closure(field, stats['mean'], stats['raw_second'])
    predicted = model_tail(model, threshold, field)
    diagnostics = {'closure': model, 'empirical_tail': stats['tail'],
                   'threshold': threshold, 'tail_comparison': '<' if field == 'm' else '>',
                   'model_tail': predicted,
                   'tail_error_model_minus_empirical': None if predicted is None else predicted - stats['tail']}
    if model['status'] == 'valid':
        distance = float(np.max(np.abs(cdf - model_cdf(model, edges))))
        if field == 'm':
            distance = max(distance, stats['atom_zero'], stats['atom_one'])
        diagnostics['max_edge_cdf_residual'] = distance
        diagnostics['model_mass_in_histogram_range'] = float(model_bins(model, edges).sum())
    return dict(stats, **diagnostics)


def save_figure(fig, out, name):
    fig.tight_layout()
    for extension in ('png', 'pdf'):
        fig.savefig(out / f'{name}.{extension}', dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_csv(path, rows, columns):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def time_plots(trace, steady_start, tau_c, out, title):
    fig, axes = plt.subplots(2, 3, figsize=(11, 5.5), sharex=True)
    for i, field in enumerate(FIELDS):
        for ax, key in zip(axes[i], ('mean', 'sigma', 'tail')):
            col = TRACE_COLUMNS.index(f'{field}_{key}')
            ax.plot(trace[:, 1], trace[:, col], lw=1)
            ax.axvline(steady_start / tau_c, color='.4', ls=':', lw=1)
            ax.set(xlabel=r'$t/\tau_c$', ylabel=f'{field}: {key}')
    fig.suptitle(title)
    save_figure(fig, out, 'stationarity_traces')


def analyse(root, out, manifest_path):
    manifest = json.loads(manifest_path.read_text())
    matches = [r for r in manifest['cases'] if root.as_posix().rstrip('/').endswith('/' + r['case'])]
    if len(matches) != 1:
        raise ValueError(f'Archive path must identify exactly one manifest case: {root}')
    row = matches[0]
    par, par_hash = read_native(root / 'parameters.json')
    verification = check_parameters(par, row)
    run_dat = manifest_path.parent / row['case'] / 'run.dat'
    if sha(run_dat) != row['run_dat_sha256']:
        raise ValueError(f'run.dat digest differs from manifest: {run_dat}')
    tau_c, sites = float(manifest['tau_c']), int(par['LX']) * int(par['LY'])
    if not np.isclose(float(par['tau_m']) / tau_c, row['tm_over_tc'], rtol=6e-6, atol=1e-12):
        raise ValueError('Manifest memory-time ratio disagrees with runtime')
    if not np.isclose(float(par['zeta_open']) / float(par['zeta']), row['activity'], rtol=6e-6, atol=1e-12):
        raise ValueError('Manifest prescribed activity disagrees with runtime')
    nsteps, ninfo, nstart = (int(par[k]) for k in ('nsteps', 'ninfo', 'nstart'))
    for key, value in (('nsteps', nsteps), ('ninfo', ninfo)):
        if value != int(row[key]):
            raise ValueError(f'Row {key} disagrees with runtime')
    if ninfo <= 0 or nsteps < nstart:
        raise ValueError('Invalid archive cadence')
    regular = [t for t in range(0, nsteps + 1, ninfo) if t >= nstart]
    expected = sorted(set(regular + [nsteps]))
    present = sorted(int(p.stem[5:]) for p in root.glob('frame*.json') if p.stem[5:].isdigit())
    if present != expected:
        raise ValueError(f'Frame set mismatch: missing={sorted(set(expected)-set(present))}, '
                         f'unexpected={sorted(set(present)-set(expected))}')
    steady_start = int(row['steady_start'])
    selected = [t for t in regular if t >= steady_start]
    if len(selected) < 4:
        raise ValueError('At least four regular steady snapshots required for two-half diagnostics')
    thresholds = {'m': float(par['mc']), 'pressure': float(par['pmem'])}
    traces, hashes = [], {}
    raw = {name: {'sum': [], 'sum2': [], 'tail_count': 0, 'atom_zero': 0,
                  'atom_one': 0, 'minimum': math.inf, 'maximum': -math.inf} for name in FIELDS}
    selected_set = set(selected)
    for t in expected:
        values, digest = native_fields(root / f'frame{t}.json', sites)
        hashes[f'frame{t}.json'] = digest
        tr = [t, t / tau_c]
        for name, x in values.items():
            tail_count = int(np.count_nonzero(x < thresholds[name] if name == 'm' else x > thresholds[name]))
            s1, s2 = float(np.sum(x)), float(np.dot(x, x))
            stats = moment_stats(s1 / sites, s2 / sites)
            tr.extend([stats[k] for k in ('mean', 'raw_second', 'variance', 'sigma')])
            tr.extend([tail_count / sites, float(x.min()), float(x.max())])
            if t in selected_set:
                acc = raw[name]
                acc['sum'].append(s1); acc['sum2'].append(s2)
                acc['tail_count'] += tail_count
                acc['atom_zero'] += int(np.count_nonzero(x == 0))
                acc['atom_one'] += int(np.count_nonzero(x == 1))
                acc['minimum'] = min(acc['minimum'], float(x.min()))
                acc['maximum'] = max(acc['maximum'], float(x.max()))
        traces.append(tr)
    trace = np.asarray(traces)
    total = sites * len(selected)
    stats = {}
    for name, acc in raw.items():
        stats[name] = dict(moment_stats(math.fsum(acc['sum']) / total, math.fsum(acc['sum2']) / total),
            tail=acc['tail_count'] / total, atom_zero=acc['atom_zero'] / total,
            atom_one=acc['atom_one'] / total, tail_count=acc['tail_count'],
            atom_zero_count=acc['atom_zero'], atom_one_count=acc['atom_one'],
            sample_count=total, minimum=acc['minimum'], maximum=acc['maximum'])
    # Shared physical bin grid gives exact summation across seeds with different ranges.
    # No pressure centering and no finite plotting bound; integer grid expands to all data.
    pressure_width = abs(float(row['expected_parameters']['zeta'])) / 256
    if pressure_width <= 0:
        raise ValueError('Nonzero baseline zeta required for common pressure bin width')
    lower = math.floor(stats['pressure']['minimum'] / pressure_width) - 1
    upper = math.ceil(stats['pressure']['maximum'] / pressure_width) + 1
    if upper - lower > 100000:
        raise ValueError('Pressure range exceeds 100000 bins; inspect unstable run before analysis')
    edges = {'pressure': np.arange(lower, upper + 1) * pressure_width,
             'm': np.unique(np.r_[np.linspace(0, 1, 513), thresholds['m']])}
    hist = {name: np.zeros(len(edges[name]) - 1, dtype=np.int64) for name in FIELDS}
    cdfs = {name: np.zeros(len(edges[name]), dtype=np.int64) for name in FIELDS}
    half_cdfs = {name: np.zeros((2, len(edges[name])), dtype=np.int64) for name in FIELDS}
    half = len(selected) // 2
    for index, t in enumerate(selected):
        values, digest = native_fields(root / f'frame{t}.json', sites)
        if digest != hashes[f'frame{t}.json']:
            raise ValueError(f'Archive changed between passes at frame{t}.json')
        for name, x in values.items():
            h = np.histogram(x, bins=edges[name])[0]
            if h.sum() != sites:
                raise ValueError(f'Histogram omitted samples from {name}')
            hist[name] += h
            c = cdf_counts(x, edges[name])
            cdfs[name] += c
            half_cdfs[name][int(index >= half)] += c
    window = trace[np.isin(trace[:, 0], selected)]
    drift, flags, arrays = {}, [], {}
    for name in FIELDS:
        cdf = cdfs[name] / total
        stats[name] = distribution_stats(stats[name], edges[name], cdf, thresholds[name], name)
        early_cdf = half_cdfs[name][0] / (half * sites)
        late_cdf = half_cdfs[name][1] / ((len(selected) - half) * sites)
        sd = max(stats[name]['sigma'], 1e-30)
        mean_col, sd_col = [TRACE_COLUMNS.index(f'{name}_{k}') for k in ('mean', 'sigma')]
        values = {'mean_shift_over_sigma': float(abs(window[:half, mean_col].mean() - window[half:, mean_col].mean()) / sd),
                  'sigma_shift_over_sigma': float(abs(window[:half, sd_col].mean() - window[half:, sd_col].mean()) / sd),
                  'early_late_edge_cdf_distance': float(np.max(abs(early_cdf - late_cdf)))}
        drift[name] = values
        for key, limit in [('mean_shift_over_sigma', .1), ('sigma_shift_over_sigma', .1),
                           ('early_late_edge_cdf_distance', .05)]:
            if values[key] > limit:
                flags.append(f'{name}: {key}={values[key]:.5g} > {limit}')
        arrays.update({f'{name}_edges': edges[name], f'{name}_count': hist[name],
                       f'{name}_density': hist[name] / total / np.diff(edges[name]),
                       f'{name}_cdf': cdf, f'{name}_model_cdf': model_cdf(stats[name]['closure'], edges[name]),
                       f'{name}_model_density': model_bins(stats[name]['closure'], edges[name]) / np.diff(edges[name]),
                       f'{name}_early_cdf': early_cdf, f'{name}_late_cdf': late_cdf})
    result = {key: row[key] for key in ('case', 'activity', 'tm_over_tc', 'seed', 'replicate')}
    result.update(schema_version=VERSION, campaign=manifest['campaign'], tau_c=tau_c,
        statistics=stats, parameters=par, parameter_verification=verification,
        source_archive=str(root.resolve()), nsteps=nsteps, ninfo=ninfo,
        frame_count=len(expected), steady_frames=len(selected), sites_per_frame=sites,
        steady_steps=selected, steady_start=steady_start, sample_count=total,
        extra_final_frame_excluded_from_statistics=nsteps not in regular,
        pressure_bin_width=pressure_width, pressure_bin_indices=[lower, upper],
        stationary=not flags, nonstationary_flags=flags, stationarity=drift,
        frame_sha256=hashes, parameters_sha256=par_hash, run_dat_sha256=row['run_dat_sha256'],
        manifest_sha256=sha(manifest_path), analysis_sha256=sha(Path(__file__)))
    np.savez_compressed(out / 'distributions.npz', **arrays)
    result['distributions_sha256'] = sha(out / 'distributions.npz')
    np.savetxt(out / 'time_series.csv', trace, delimiter=',', header=','.join(TRACE_COLUMNS), comments='')
    result['time_series_sha256'] = sha(out / 'time_series.csv')
    dump(out / 'distribution_result.json', result)
    time_plots(trace, steady_start, tau_c, out,
               f"a={row['activity']:g}, g={row['tm_over_tc']:g}, seed={row['seed']}"
               + (' — drift flagged' if flags else ''))
    print(json.dumps({'case': row['case'], 'steady_frames': len(selected),
                      'stationary': not flags, 'statistics': stats}), flush=True)


def pad_pressure(array, left, right, is_cdf=False):
    return np.pad(array, (left, right), mode='constant', constant_values=(0, 1) if is_cdf else 0)


def summarize_group(rows, root, out):
    if len(rows) != 3 or len({r['seed'] for r in rows}) != 3:
        raise ValueError('Every (activity, tau_m) group must contain three distinct seeds')
    rows = sorted(rows, key=lambda r: r['seed'])
    widths = [r['pressure_bin_width'] for r in rows]
    if not np.allclose(widths, widths[0], rtol=0, atol=0):
        raise ValueError('Incompatible physical pressure bin grids')
    lower = min(r['pressure_bin_indices'][0] for r in rows)
    upper = max(r['pressure_bin_indices'][1] for r in rows)
    result = {'activity': rows[0]['activity'], 'tm_over_tc': rows[0]['tm_over_tc'],
              'seeds': [r['seed'] for r in rows], 'cases': [r['case'] for r in rows],
              'stationary': all(r['stationary'] for r in rows),
              'nonstationary_cases': [r['case'] for r in rows if not r['stationary']],
              'seed_weight': 1 / 3, 'statistics': {}}
    arrays = {}
    distributions = []
    for row in rows:
        path = root / row['case'] / 'distributions.npz'
        if sha(path) != row['distributions_sha256']:
            raise ValueError(f'Distribution digest mismatch: {path}')
        with np.load(path) as data:
            distributions.append({key: data[key] for key in data.files})
    for name in FIELDS:
        edges = distributions[0]['m_edges'] if name == 'm' else np.arange(lower, upper + 1) * widths[0]
        seed_density, seed_cdf = [], []
        for row, data in zip(rows, distributions):
            if name == 'm':
                if not np.array_equal(data['m_edges'], edges):
                    raise ValueError('Incompatible memory histogram grids')
                density, cdf = data['m_density'], data['m_cdf']
            else:
                left, right = row['pressure_bin_indices'][0] - lower, upper - row['pressure_bin_indices'][1]
                density = pad_pressure(data['pressure_density'], left, right)
                cdf = pad_pressure(data['pressure_cdf'], left, right, is_cdf=True)
            seed_density.append(density); seed_cdf.append(cdf)
        density, cdf = np.mean(seed_density, axis=0), np.mean(seed_cdf, axis=0)
        if not np.isclose(np.sum(density * np.diff(edges)), 1, atol=1e-12):
            raise ValueError('Equal-seed pooled PDF normalization failed')
        values = [r['statistics'][name] for r in rows]
        stats = moment_stats(np.mean([s['mean'] for s in values]), np.mean([s['raw_second'] for s in values]))
        for metric in ('tail', 'atom_zero', 'atom_one'):
            stats[metric] = float(np.mean([s[metric] for s in values]))
        stats['minimum'] = min(s['minimum'] for s in values)
        stats['maximum'] = max(s['maximum'] for s in values)
        threshold = values[0]['threshold']
        if any(s['threshold'] != threshold for s in values):
            raise ValueError('Inconsistent thresholds within a group')
        stats = distribution_stats(stats, edges, cdf, threshold, name)
        stats['seed_sd'] = {key: float(np.std([s[key] for s in values], ddof=1))
                            for key in ('mean', 'raw_second', 'variance', 'sigma', 'tail')}
        stats['mean_seed_sigma'] = float(np.mean([s['sigma'] for s in values]))
        result['statistics'][name] = stats
        arrays.update({f'{name}_edges': edges, f'{name}_density': density,
            f'{name}_density_seed_sd': np.std(seed_density, axis=0, ddof=1),
            f'{name}_cdf': cdf, f'{name}_cdf_seed_sd': np.std(seed_cdf, axis=0, ddof=1),
            f'{name}_seed_density': np.asarray(seed_density), f'{name}_seed_cdf': np.asarray(seed_cdf),
            f'{name}_model_density': model_bins(stats['closure'], edges) / np.diff(edges),
            f'{name}_model_cdf': model_cdf(stats['closure'], edges)})
    stem = f"g{result['tm_over_tc']:g}_a{result['activity']:g}".replace('.', 'p')
    np.savez_compressed(out / f'{stem}.npz', seeds=result['seeds'], **arrays)
    result['distribution_file'] = f'{stem}.npz'
    return result, arrays


def group_plots(groups, arrays, out):
    times = sorted({g['tm_over_tc'] for g in groups})
    for time in times:
        subset = sorted([g for g in groups if g['tm_over_tc'] == time], key=lambda g: g['activity'])
        for name in FIELDS:
            fig, axes = plt.subplots(3, len(subset), figsize=(3.4 * len(subset), 8), squeeze=False)
            for column, row in enumerate(subset):
                values = arrays[(time, row['activity'])]
                stats = row['statistics'][name]
                edges = values[f'{name}_edges']
                centers = (edges[:-1] + edges[1:]) / 2
                density, spread = values[f'{name}_density'], values[f'{name}_density_seed_sd']
                axes[0, column].fill_between(centers, np.maximum(0, density - spread), density + spread,
                                             color='#2166ac', alpha=.18, label='seed SD')
                axes[0, column].stairs(density, edges, color='#2166ac', lw=1.1, label='native data')
                axes[0, column].stairs(values[f'{name}_model_density'], edges,
                                       color='#b35806', ls='--', lw=1.2, label=stats['closure']['family'])
                for seed_cdf in values[f'{name}_seed_cdf']:
                    axes[1, column].plot(edges, seed_cdf, color='#2166ac', alpha=.25, lw=.6)
                    axes[2, column].plot(edges, seed_cdf - values[f'{name}_model_cdf'], color='#2166ac', alpha=.25, lw=.6)
                axes[1, column].plot(edges, values[f'{name}_cdf'], color='#2166ac')
                axes[1, column].plot(edges, values[f'{name}_model_cdf'], color='#b35806', ls='--')
                axes[2, column].plot(edges, values[f'{name}_cdf'] - values[f'{name}_model_cdf'], color='#2166ac')
                axes[2, column].axhline(0, color='.5', lw=.6)
                for ax in axes[:, column]:
                    ax.axvline(stats['threshold'], color='.4', ls=':', lw=.8)
                    ax.set_xlim(edges[0], edges[-1])
                    ax.set_xlabel('$m$' if name == 'm' else '$P$')
                axes[0, column].set_title(f"a={row['activity']:g}" + ('; drift flagged' if not row['stationary'] else ''))
                if stats['closure']['status'] != 'valid':
                    axes[0, column].text(.05, .85, 'Continuous model undefined', transform=axes[0, column].transAxes, fontsize=8)
            for ax, label in zip(axes[:, 0], ('Probability density', r'$\Pr(X\leq x)$', 'Data CDF − model CDF')):
                ax.set_ylabel(label)
            axes[0, 0].legend(fontsize=8)
            fig.suptitle(f"{'Memory / Beta' if name == 'm' else 'Pressure / Gaussian'}; "
                         + rf'$\tau_m/\tau_c={time:g}$' + '; equal mean and second moment')
            save_figure(fig, out, f'{name}_g{time:g}'.replace('.', 'p'))


def moment_plots(groups, out):
    for horizontal in ('tm_over_tc', 'activity'):
        other = 'activity' if horizontal == 'tm_over_tc' else 'tm_over_tc'
        fig, axes = plt.subplots(2, 3, figsize=(12, 6))
        for value in sorted({g[other] for g in groups}):
            subset = sorted([g for g in groups if g[other] == value], key=lambda g: g[horizontal])
            x = [g[horizontal] for g in subset]
            for i, name in enumerate(FIELDS):
                for ax, metric in zip(axes[i], ('mean', 'raw_second', 'tail')):
                    y = [g['statistics'][name][metric] for g in subset]
                    sd = [g['statistics'][name]['seed_sd'][metric] for g in subset]
                    label = f"{'a' if other == 'activity' else 'g'}={value:g}"
                    ax.errorbar(x, y, yerr=sd, marker='o', ms=3, lw=1, capsize=2, label=label)
                    for xx, yy, g in zip(x, y, subset):
                        if not g['stationary']:
                            ax.plot(xx, yy, marker='x', color='black', ms=7, linestyle='none')
                    ax.set(ylabel=f'{name}: {metric}', xlabel=r'$\tau_m/\tau_c$' if horizontal == 'tm_over_tc' else '$a$')
                    if horizontal == 'tm_over_tc':
                        ax.set_xscale('log')
        axes[0, 0].legend(fontsize=8, ncol=2)
        fig.suptitle('Equal-seed mean ± seed SD; × = stationarity diagnostic flagged')
        save_figure(fig, out, f'moments_vs_{horizontal}')


def summarize(root, out, manifest_path):
    manifest = json.loads(manifest_path.read_text())
    manifest_hash = sha(manifest_path)
    rows, missing, invalid = [], [], []
    for case in manifest['cases']:
        path = root / case['case'] / 'distribution_result.json'
        if not path.exists():
            missing.append(case['case']); continue
        try:
            row = json.loads(path.read_text())
            if row['manifest_sha256'] != manifest_hash or row['run_dat_sha256'] != case['run_dat_sha256']:
                raise ValueError('Manifest or input provenance mismatch')
            if row['analysis_sha256'] != sha(Path(__file__)):
                raise ValueError('Analysis code changed since per-case extraction')
            if any(row[key] != case[key] for key in ('case', 'activity', 'tm_over_tc', 'seed', 'replicate')):
                raise ValueError('Result identity differs from manifest')
            rows.append(row)
        except (ValueError, KeyError) as error:
            invalid.append({'case': case['case'], 'error': str(error)})
    report = {'schema_version': VERSION, 'campaign': manifest['campaign'],
              'expected': len(manifest['cases']), 'available': len(rows),
              'missing': missing, 'invalid': invalid,
              'nonstationary': [r['case'] for r in rows if not r['stationary']],
              'manifest_sha256': manifest_hash, 'analysis_sha256': sha(Path(__file__)),
              'cases': rows, 'groups': []}
    dump(out / 'summary.json', report)
    if missing or invalid:
        raise RuntimeError(f'Missing {len(missing)} / invalid {len(invalid)} cases; summary.json retains report')
    grouped = {}
    for row in rows:
        grouped.setdefault((row['tm_over_tc'], row['activity']), []).append(row)
    groups, arrays = [], {}
    for key, members in sorted(grouped.items()):
        group, data = summarize_group(members, root, out)
        groups.append(group); arrays[key] = data
    report['groups'] = groups
    dump(out / 'summary.json', report)
    per_seed, per_group = [], []
    for source, destination in ((rows, per_seed), (groups, per_group)):
        for row in source:
            for name, stats in row['statistics'].items():
                flat = {key: row.get(key, '') for key in ('case', 'activity', 'tm_over_tc', 'seed', 'replicate', 'stationary')}
                flat.update(field=name, **{k: v for k, v in stats.items() if not isinstance(v, dict)})
                flat.update({f'seed_sd_{key}': val for key, val in stats.get('seed_sd', {}).items()})
                flat.update(model_status=stats['closure']['status'], model=stats['closure']['family'],
                            alpha=stats['closure'].get('alpha', ''), beta=stats['closure'].get('beta', ''))
                destination.append(flat)
    for name, values in (('per_seed_statistics', per_seed), ('group_statistics', per_group)):
        columns = list(dict.fromkeys(key for row in values for key in row))
        write_csv(out / f'{name}.csv', values, columns)
    group_plots(groups, arrays, out)
    moment_plots(groups, out)
    (out / 'README.md').write_text('''# Poster distribution calibration

Native fields are `pressure` (physical mechanical pressure in the simulation gauge)
and `m`, read directly from JSON. No centering, clipping, smoothing, optimized fit,
or video quantization. This is an analysis summary, not final poster styling.

Each seed averages every grid point and every regular snapshot at/after its declared
`steady_start`; an off-cadence final frame is checked but excluded from statistics.
The first moment and raw second moment are calculated directly from native numbers,
not from histograms. Variance is raw second moment minus mean squared. Gaussian
pressure and Beta memory parameters are fixed by those moments. Analytic model
moments and residuals are exported so matching can be verified independently.
Degenerate/Bernoulli-limit moments are reported without inventing a continuous PDF.

The three seeds have equal weight even when their snapshot counts differ. Group
moments are the mean of seed raw moments; group variance includes between-seed mean
variation. Every seed is retained, including stationarity flags. Error bars and bands
show sample SD across seeds, not confidence intervals; lattice sites and consecutive
frames are not treated as independent replicates. `sigma` for a group means the
pooled distribution sigma, while `mean_seed_sigma` means the average of seed sigmas.

Memory bins span [0,1], with m_c inserted as an edge. Pressure bins use a shared
physical width zeta/256 and extend beyond both native extrema, preserving all tails.
Bin-integrated model densities avoid Beta endpoint singularities. Gaussian tails
outside the displayed native range are neither renormalized nor hidden in the
export: `model_mass_in_histogram_range` records the captured model mass.

CDFs are counted directly from fields with <= at every edge. Threshold statistics
are counted separately with the actual strict inequalities m < m_c and P > p_c.
Atoms at m=0 and m=1 are retained. CDF residuals are edge diagnostics, not independent-
site KS tests or p-values. NPZ files retain every seed PDF/CDF, seed SD, equal-seed
PDF/CDF and model curves; JSON/CSV retain native moments and exact threshold counts.

Two-half stationarity diagnostics flag a mean shift >0.1 pooled sigma, a snapshot
sigma shift >0.1 pooled sigma, or an edge-CDF change >0.05. Passing these checks is
not proof of stationarity. Flags never remove a seed. Each case has full-time traces
in `time_series.csv` and `stationarity_traces.png/pdf`, plus runtime parameter, input,
frame and analysis hashes. Missing/invalid cases cause a nonzero exit after saving
a report; incomplete group figures are not presented as complete.
''')
    print(json.dumps({'expected': len(rows), 'groups': len(groups),
                      'nonstationary': report['nonstationary']}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument('--summary', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(STYLE)
    try:
        (summarize if args.summary else analyse)(args.root, args.out, args.manifest)
    except Exception as error:
        dump(args.out / 'analysis_failure.json', {'error': type(error).__name__, 'message': str(error),
             'root': str(args.root), 'manifest': str(args.manifest), 'summary_mode': args.summary,
             'analysis_sha256': sha(Path(__file__))})
        raise
    failure = args.out / 'analysis_failure.json'
    if failure.exists():
        failure.unlink()


if __name__ == '__main__':
    main()
