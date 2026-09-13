#!/usr/bin/env python3
"""Pressure-threshold pulse diagnostics and conservative exponential relaxation fits.

Per case: INPUT_RAW OUTPUT_ANALYSIS -> pulse_case.json, response_series.npz.
Campaign: INPUT_RESULTS OUTPUT_SUMMARY --summary --manifest manifest.json.

Only full-domain C++ response.csv statistics are used. A control and its pulse
are aligned on identical absolute simulation steps, without interpolation.
Stationarity checks are finite-window diagnostics, not stationarity proofs.
Fits have exactly two parameters, A and tau, and no fitted constant offset.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re

import numpy as np

TC = 674.3290333006435
FIELDS = ('step', 'pmem_effective', 'pulse_on', 'chi_mean', 'chi_std',
          'm_mean', 'm_std', 'u_rms', 'P_mean', 'P_std', 'source_mean',
          'source_base_mean', 'chi_target_mean')
CRITERIA = {
    'baseline_windows_tc': {'far_from_transition': 50., 'near_transition': 200.},
    'near_transition_tm_over_tc': [7., 9.],
    'terminal_window_max_tc': 100.,
    'drift_limits': {'chi_mean': {'half_change': .015, 'block_range': .03},
                     'm_mean': {'half_change': .008, 'block_range': .016},
                     'u_rms': {'half_change_relative': .10, 'block_range_relative': .20}},
    'control_baseline_to_tail_absolute_limits': {'chi_mean': .03, 'm_mean': .015},
    'return_absolute_tolerances': {'chi_mean': .015, 'm_mean': .008},
    'return_noise_multiplier': 2.,
    'switch_minimum_chi_separation': .20,
    'fit_bin_tc': .1,
    'minimum_peak_noise_ratio': 5.,
    'minimum_fit_amplitude_noise_ratio': 4.,
    'minimum_fit_bins': 12,
    'minimum_fit_r_squared': .80,
    'minimum_window_in_fitted_tau': 1.5,
    'maximum_relative_start_sensitivity': .25,
}


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def parameters(root):
    data = json.loads((root / 'parameters.json').read_text())['data']
    return {key: item['value'] for key, item in data.items()}


def blocks(t, values, start, end, n=5):
    edges = np.linspace(start, end, n + 1)
    return np.array([np.mean(values[(t >= a) & (t < b)])
                     if np.any((t >= a) & (t < b)) else np.nan
                     for a, b in zip(edges[:-1], edges[1:])])


def diagnostics(t, series, start, end):
    """Unweighted means of dense fixed-cadence samples; event rows are negligible."""
    result = {'start_tc': float(start), 'end_tc': float(end),
              'duration_tc': float(end - start), 'metrics': {}, 'flags': []}
    if end <= start:
        result['flags'].append('empty_window')
        result['usable'] = False
        return result
    for name in ('chi_mean', 'm_mean', 'u_rms'):
        b = blocks(t, series[name], start, end)
        h = blocks(t, series[name], start, end, 2)
        if not np.isfinite(np.r_[b, h]).all():
            result['flags'].append(name + '_insufficient_samples')
            continue
        avg, drift, span = float(np.mean(b)), float(h[1] - h[0]), float(np.ptp(b))
        relative = name == 'u_rms'
        scale = max(abs(avg), 1e-12) if relative else 1.
        limits = CRITERIA['drift_limits'][name]
        dlim = limits['half_change_relative' if relative else 'half_change']
        rlim = limits['block_range_relative' if relative else 'block_range']
        result['metrics'][name] = {'mean': avg, 'blocks': b.tolist(),
                                  'block_std': float(np.std(b)),
                                  'half_change': drift, 'block_range': span,
                                  'relative_half_change': drift / scale if relative else None}
        if abs(drift) / scale > dlim:
            result['flags'].append(name + '_drift')
        if span / scale > rlim:
            result['flags'].append(name + '_block_variation')
        if relative and avg <= 1e-8:
            result['flags'].append('flow_not_developed')
    result['usable'] = not result['flags']
    return result


def reduce_case(root, out):
    p = parameters(root)
    lines = '\n'.join(s for s in (root / 'response.csv').read_text().splitlines()
                      if not s.startswith('#'))
    data = np.atleast_1d(np.genfromtxt(io.StringIO(lines), delimiter=',', names=True))
    if data.dtype.names is None or not set(FIELDS).issubset(data.dtype.names):
        raise ValueError('response.csv is missing required columns')
    if len(data) < 10 or any(not np.isfinite(data[name]).all() for name in FIELDS):
        raise ValueError('Insufficient or non-finite response records')
    steps = data['step']
    if np.any(steps != np.rint(steps)) or np.any(np.diff(steps) <= 0):
        raise ValueError('Response steps must be integral and strictly increasing')
    nsteps, freeze = int(p['nsteps']), int(p.get('chi_freeze_steps', 0))
    start, duration = int(p['pmem_pulse_start']), int(p['pmem_pulse_steps'])
    cadence = int(p['nresponse'])
    if cadence <= 0 or not (0 <= freeze < start < nsteps) or duration < 0:
        raise ValueError('Invalid pulse/feedback/response clock parameters')
    if start + duration > nsteps:
        raise ValueError('Pulse extends beyond the observation window')
    expected = set(range(0, nsteps + 1, cadence)) | {freeze, start, nsteps}
    if duration > 0:
        expected.add(start + duration)
    if set(map(int, steps)) != expected:
        missing = sorted(expected - set(map(int, steps)))
        extra = sorted(set(map(int, steps)) - expected)
        raise ValueError(f'Incomplete/incorrect response clock: missing={missing[:8]}, extra={extra[:8]}')
    on = (duration > 0) & (steps >= start) & (steps < start + duration)
    if not np.array_equal(data['pulse_on'], on.astype(float)):
        raise ValueError('pulse_on disagrees with the half-open pulse interval')
    expected_pmem = np.where(on, float(p['pmem_pulse_value']), float(p['pmem']))
    if not np.allclose(data['pmem_effective'], expected_pmem, atol=1e-9, rtol=5e-6):
        raise ValueError('Effective pmem differs from requested baseline/pulse values')
    for name in ('chi_mean', 'm_mean', 'source_mean', 'source_base_mean', 'chi_target_mean'):
        if np.min(data[name]) < -1e-6 or np.max(data[name]) > 1 + 1e-6:
            raise ValueError(name + ' outside [0,1]')
    for name in ('chi_std', 'm_std', 'u_rms', 'P_std'):
        if np.min(data[name]) < -1e-12:
            raise ValueError(name + ' cannot be negative')
    match = re.fullmatch(r'tm([0-9]+(?:p[0-9]+)?)_init(patch|0|1)_rep([1-4])', root.name)
    if match is None:
        raise ValueError('Missing tau_m, initialization or replicate in case name')
    tm, initialization, replicate = (float(match[1].replace('p', '.')), match[2], int(match[3]))
    if not np.isclose(float(p['tau_m']) / TC, tm, rtol=5e-6):
        raise ValueError('Case name and serialized tau_m disagree')
    clock = steps / TC
    near_low, near_high = CRITERIA['near_transition_tm_over_tc']
    region = 'near_transition' if near_low <= tm <= near_high else 'far_from_transition'
    baseline_width = CRITERIA['baseline_windows_tc'][region]
    baseline_start = max(freeze / TC, start / TC - baseline_width)
    baseline = diagnostics(clock, data, baseline_start, start / TC)
    baseline['requested_window_tc'] = baseline_width
    if baseline['duration_tc'] < baseline_width - .05:
        baseline['flags'].append('baseline_window_too_short')
        baseline['usable'] = False
    full_feedback = diagnostics(clock, data, freeze / TC, start / TC)
    full_feedback['interpretation'] = (
        'Descriptive evolution from feedback release to pulse onset, including initial transients; '
        'these flags do not gate baseline_usable or control_usable.')
    full_feedback['gates_baseline_usable'] = False
    off = start + duration
    tail_width = min(CRITERIA['terminal_window_max_tc'], (nsteps - off) / TC / 2.)
    terminal = diagnostics(clock, data, nsteps / TC - tail_width,
                           np.nextafter(nsteps / TC, np.inf))
    control_flags = []
    if duration == 0:
        control_flags.extend(terminal['flags'])
        for name, limit in CRITERIA['control_baseline_to_tail_absolute_limits'].items():
            if name not in baseline['metrics'] or name not in terminal['metrics']:
                continue
            shift = terminal['metrics'][name]['mean'] - baseline['metrics'][name]['mean']
            if abs(shift) > limit:
                control_flags.append(name + '_baseline_to_tail_drift')
        if all('u_rms' in d['metrics'] for d in (baseline, terminal)):
            b = baseline['metrics']['u_rms']['mean']
            if abs(terminal['metrics']['u_rms']['mean'] - b) / max(abs(b), 1e-12) > .20:
                control_flags.append('u_rms_baseline_to_tail_drift')
    row = {'case_name': root.name, 'kind': 'pulse' if duration else 'control',
           'tm_over_tc': tm, 'initialization': initialization, 'replicate': replicate,
           'seed': int(p['seed']), 'L': int(p['LX']), 'mc': float(p['mc']),
           'pmem': float(p['pmem']), 'pmem_pulse_value': float(p['pmem_pulse_value']),
           'tau_chi': float(p['tau_chi']), 'tau_c_reference': TC,
           'nsteps': nsteps, 'preparation_steps': freeze,
           'pulse_start_steps': start, 'pulse_duration_steps': duration,
           'pulse_end_steps': off, 'nresponse': cadence, 'records': len(data),
           'complete': True, 'baseline': baseline, 'full_feedback_window': full_feedback,
           'terminal': terminal,
           'control_post_flags': sorted(set(control_flags)),
           'baseline_usable': baseline['usable'],
           'control_usable': baseline['usable'] and not control_flags if duration == 0 else None,
           'criteria': CRITERIA,
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'interpretation': 'Finite-window screening only; correlated samples are not independent trials.'}
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / 'response_series.npz', **{name: data[name] for name in FIELDS})
    write_json(out / 'pulse_case.json', row)
    print(json.dumps({'case': root.name, 'complete': True, 'baseline_usable': row['baseline_usable'],
                      'control_post_flags': row['control_post_flags']}))
    return row


def bin_series(t, y, width=.1):
    indices = np.floor(np.maximum(t, 0.) / width + 1e-10).astype(int)
    unique, inverse = np.unique(indices, return_inverse=True)
    count = np.bincount(inverse)
    return (np.bincount(inverse, weights=t) / count,
            np.bincount(inverse, weights=y) / count)


def exponential_fit(t, response):
    """Fit signed late tails and reject weak/unresolved/multimode windows.

    No time-sample-derived confidence interval is reported. Independent seeds,
    and sensitivity to the fit start, provide the reported spread instead.
    """
    from scipy.optimize import curve_fit
    t, y = bin_series(np.asarray(t), np.asarray(response), CRITERIA['fit_bin_tc'])
    report = {'status': 'insufficient_samples', 'tau_tc': None, 'candidates': []}
    if len(t) < 30:
        return report
    tail = y[max(len(y) * 3 // 4, 0):]
    noise = max(1e-8, float(1.4826 * np.median(np.abs(tail - np.median(tail)))))
    # Local averaging is used only to locate the end of the initial excursion.
    # The fit itself uses unsmoothed bin averages.
    smooth = np.convolve(y, np.ones(3) / 3., mode='same')
    smooth[0], smooth[-1] = y[0], y[-1]
    peak = int(np.argmax(np.abs(smooth)))
    sign = float(np.sign(smooth[peak]))
    amplitude = float(abs(smooth[peak]))
    report.update({'noise_floor': noise, 'peak_response': sign * amplitude,
                   'peak_time_tc': float(t[peak]), 'peak_noise_ratio': amplitude / noise})
    if amplitude < CRITERIA['minimum_peak_noise_ratio'] * noise:
        report['status'] = 'signal_too_weak'
        return report
    positive = sign * y
    envelope = sign * smooth
    after = np.flatnonzero((np.arange(len(t)) > peak) & (envelope < amplitude / np.e))
    if not len(after):
        report['status'] = 'decay_not_resolved'
        return report
    rough_tau = max(t[after[0]] - t[peak], .1)
    # End the window on sustained arrival at the noise floor; isolated crossings
    # do not truncate an otherwise resolved tail.
    hold = max(3, int(np.ceil(min(2., max(.3, rough_tau / 4.)) / .1)))
    endpoint = len(t) - 1
    for j in range(peak + 1, len(t) - hold):
        if np.all(envelope[j:j + hold] < 2 * noise):
            endpoint = j
            break
    fit_end = float(t[endpoint])
    report.update({'rough_tau_tc': float(rough_tau), 'fit_end_tc': fit_end})

    def model(tt, aa, tau):
        return aa * np.exp(-tt / tau)

    for offset in (0., .25, .5, 1., 1.5, 2., 3., 4.):
        requested_start = t[peak] + offset * rough_tau
        mask = (t >= requested_start) & (t <= fit_end)
        if mask.sum() < CRITERIA['minimum_fit_bins']:
            continue
        tt, yy = t[mask], positive[mask]
        tt0 = tt - tt[0]
        candidate = {'start_tc': float(tt[0]), 'end_tc': float(tt[-1]),
                     'bins': int(mask.sum()), 'accepted': False}
        try:
            best, _ = curve_fit(model, tt0, yy,
                                p0=(max(float(yy[0]), noise), rough_tau),
                                bounds=([0., .005], [np.inf, max(1., 10 * (t[-1] - t[0]))]),
                                maxfev=10000)
            amp, tau = map(float, best)
            residual = yy - model(tt0, *best)
            sse = float(residual @ residual)
            sst = float(np.sum((yy - np.mean(yy)) ** 2))
            r2 = 1 - sse / sst if sst > 1e-24 else -1.
            candidate.update({'tau_tc': tau, 'amplitude_at_start': sign * amp,
                              'r_squared': r2, 'rmse': float(np.sqrt(sse / len(yy))),
                              'window_in_tau': float(tt0[-1] / tau)})
            candidate['accepted'] = bool(
                amp >= CRITERIA['minimum_fit_amplitude_noise_ratio'] * noise
                and r2 >= CRITERIA['minimum_fit_r_squared']
                and tt0[-1] / tau >= CRITERIA['minimum_window_in_fitted_tau'])
        except (RuntimeError, ValueError, FloatingPointError) as exc:
            candidate['error'] = str(exc)
        report['candidates'].append(candidate)
    good = [r for r in report['candidates'] if r['accepted']]
    report['status'] = 'no_resolved_single_exponential_window'
    if len(good) < 2:
        return report
    # Use the latest pair with stable tau: fast early modes can have decayed by
    # this point. Retain every candidate so the choice is inspectable.
    for left, right in reversed(list(zip(good[:-1], good[1:]))):
        spread = abs(left['tau_tc'] - right['tau_tc']) / np.mean([left['tau_tc'], right['tau_tc']])
        if spread <= CRITERIA['maximum_relative_start_sensitivity']:
            report.update({'status': 'resolved', 'tau_tc': right['tau_tc'],
                           'fit_start_tc': right['start_tc'], 'fit_end_tc': right['end_tc'],
                           'amplitude_at_start': right['amplitude_at_start'],
                           'r_squared': right['r_squared'],
                           'start_sensitivity_relative': float(spread),
                           'sensitivity_tau_range_tc': [min(left['tau_tc'], right['tau_tc']),
                                                        max(left['tau_tc'], right['tau_tc'])]})
            return report
    report['status'] = 'fit_start_sensitive'
    return report


def pair_response(pulse, control, pseries, cseries, alternate_control=None, allow_drift=False):
    result = {'classification': 'invalid_pair', 'fits': {}, 'flags': []}
    for key in ('tm_over_tc', 'initialization', 'replicate', 'seed', 'L', 'mc', 'pmem',
                'tau_chi', 'preparation_steps', 'pulse_start_steps'):
        if pulse[key] != control[key]:
            result['flags'].append('mismatched_' + key)
    if control['kind'] != 'control' or pulse['kind'] != 'pulse':
        result['flags'].append('wrong_pair_kind')
    if result['flags']:
        return result, None
    absolute, pi, ci = np.intersect1d(pseries['step'], cseries['step'], return_indices=True)
    before = absolute < pulse['pulse_start_steps']
    if before.sum() < 10:
        result['flags'].append('missing_paired_history')
    # A deterministic, identical prehistory is a stronger check than equal seed.
    history_fields = ('chi_mean', 'chi_std', 'm_mean', 'm_std', 'u_rms', 'P_mean', 'P_std')
    differences = {k: float(np.max(np.abs(pseries[k][pi[before]] - cseries[k][ci[before]])))
                   for k in history_fields} if before.any() else {}
    result['prehistory_max_absolute_differences'] = differences
    if any(v > 1e-9 for v in differences.values()):
        result['flags'].append('prehistory_not_identical')
    if result['flags']:
        return result, None
    after = absolute >= pulse['pulse_end_steps']
    if after.sum() < 30 or absolute[-1] < pulse['nsteps']:
        result['flags'].append('control_does_not_cover_pulse_recovery')
        return result, None
    t = (absolute[after] - pulse['pulse_end_steps']) / TC
    series = {'time_after_pulse_tc': t, 'step': absolute[after]}
    for name in ('chi_mean', 'm_mean'):
        series['delta_' + name] = pseries[name][pi[after]] - cseries[name][ci[after]]
    result['baseline_stationary_by_legacy_checks'] = bool(pulse['baseline_usable'] and control['control_usable'])
    result['response_interpretation'] = ('Fixed-age paired response; drift diagnostics do not veto measurement. '
        'A fitted decay time alone does not establish stationary-state stability.' if allow_drift else
        'Stationarity-screened return to the original state.')
    if (not pulse['baseline_usable'] or not control['control_usable']) and not allow_drift:
        result['classification'] = 'baseline_drift'
        result['flags'].extend(['pulse_baseline:' + f for f in pulse['baseline']['flags']])
        result['flags'].extend(['control_baseline:' + f for f in control['baseline']['flags']])
        result['flags'].extend(['control_post:' + f for f in control['control_post_flags']])
        return result, series
    tail = t >= t[-1] - min(CRITERIA['terminal_window_max_tc'], t[-1] / 2.)
    shifts, tolerances = {}, {}
    for name in ('chi_mean', 'm_mean'):
        shifts[name] = float(np.mean(series['delta_' + name][tail]))
        base_tolerance = CRITERIA['return_absolute_tolerances'][name]
        noise = control['terminal']['metrics'][name]['block_std']
        tolerances[name] = max(base_tolerance, CRITERIA['return_noise_multiplier'] * noise)
    result.update({'terminal_paired_shift': shifts, 'return_tolerances': tolerances})
    same_state = all(abs(shifts[k]) <= tolerances[k] for k in shifts)
    if allow_drift:
        result['diagnostic_flags'] = (['pulse_baseline:' + f for f in pulse['baseline']['flags']]
            + ['control_baseline:' + f for f in control['baseline']['flags']]
            + ['control_post:' + f for f in control['control_post_flags']]
            + ['pulse_terminal:' + f for f in pulse['terminal']['flags']])
        result['classification'] = 'paired_response_decayed' if same_state else 'persistent_paired_difference'
        result['control_terminal_chi'] = control['terminal']['metrics']['chi_mean']['mean']
        result['pulse_terminal_chi'] = pulse['terminal']['metrics']['chi_mean']['mean']
        if same_state:
            for name in ('chi_mean', 'm_mean'):
                result['fits'][name] = exponential_fit(t, series['delta_' + name])
            if all(f['status'] == 'signal_too_weak' for f in result['fits'].values()):
                result['classification'] = 'signal_too_weak'
        return result, series
    if same_state and pulse['terminal']['usable']:
        result['classification'] = 'returned_to_original_state'
        for name in ('chi_mean', 'm_mean'):
            result['fits'][name] = exponential_fit(t, series['delta_' + name])
        if all(f['status'] == 'signal_too_weak' for f in result['fits'].values()):
            result['classification'] = 'signal_too_weak'
        return result, series
    # A switch requires a stable terminal window AND evidence of a distinct
    # branch. Prefer the opposite-initialization control with this same seed.
    original = control['terminal']['metrics']['chi_mean']['mean']
    final = pulse['terminal']['metrics']['chi_mean']['mean']
    minimum_gap = CRITERIA['switch_minimum_chi_separation']
    target = None
    if alternate_control and alternate_control['control_usable']:
        target = alternate_control['terminal']['metrics']['chi_mean']['mean']
        result['alternate_control_chi'] = target
    matches_other = (target is not None and abs(target - original) >= minimum_gap
                     and abs(final - target) <= max(.05, 2 * tolerances['chi_mean']))
    opposite_extreme = ((original >= .65 and final <= .35)
                        or (original <= .35 and final >= .65))
    if pulse['terminal']['usable'] and abs(final - original) >= minimum_gap and (matches_other or opposite_extreme):
        result['classification'] = 'switched_branch' if matches_other else 'switch_candidate'
        # This is a coarse first persistent entry time, not an escape-rate fit.
        midpoint = (original + (target if target is not None else final)) / 2.
        direction = np.sign(final - original)
        from_onset = absolute >= pulse['pulse_start_steps']
        onset_clock = (absolute[from_onset] - pulse['pulse_start_steps']) / TC
        bt, bx = bin_series(onset_clock, pseries['chi_mean'][pi[from_onset]], 2.)
        persistent = direction * (bx - midpoint) > 0
        dwell = max(5, int(np.ceil(pulse['tm_over_tc'] / 2.)))
        entry = next((j for j in range(max(0, len(bt) - dwell + 1))
                      if np.all(persistent[j:j + dwell])), None)
        duration_tc = pulse['pulse_duration_steps'] / TC
        result['switch_entry_since_pulse_onset_tc'] = float(bt[entry]) if entry is not None else None
        result['switch_entry_after_pulse_tc'] = float(bt[entry] - duration_tc) if entry is not None else None
        result['switch_entry_resolution_tc'] = 2.
        result['switch_entry_dwell_tc'] = 2. * dwell
    else:
        result['classification'] = 'not_recovered'
    result['flags'].extend(['pulse_terminal:' + f for f in pulse['terminal']['flags']])
    return result, series


def locate_case(root, name):
    for parent in (root / name / 'analysis', root / name):
        if (parent / 'pulse_case.json').is_file():
            return parent
    return None


def aggregate(root, out, manifest_path, make_figures=True, allow_drift=False):
    manifest = json.loads(manifest_path.read_text())
    rows, directories, missing, invalid = {}, {}, [], []
    fixed = manifest['fixed']
    names = [item['case'] for item in manifest['cases']]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate case identifiers in manifest')
    for item in manifest['cases']:
        case = item['case']
        directory = locate_case(root, case)
        if directory is None:
            missing.append(case)
            continue
        try:
            row = json.loads((directory / 'pulse_case.json').read_text())
            if row['case_name'] != Path(case).name:
                raise ValueError('Wrong case identity')
            for key in ('initialization', 'replicate', 'kind', 'pulse_start_steps',
                        'pulse_duration_steps', 'nsteps', 'preparation_steps'):
                if str(row[key]) != str(item[key]):
                    raise ValueError('Wrong ' + key)
            if not np.isclose(row['tm_over_tc'], item['tm_over_tc'], rtol=1e-10):
                raise ValueError('Wrong tau_m coordinate')
            for key in ('mc', 'pmem', 'pmem_pulse_value', 'tau_chi', 'L'):
                if not np.isclose(row[key], fixed[key], rtol=5e-6, atol=1e-9):
                    raise ValueError('Wrong ' + key)
            if not np.isclose(row['tau_c_reference'], fixed['tau_c'], rtol=1e-10):
                raise ValueError('Wrong reference tau_c')
            if not (directory / 'response_series.npz').is_file():
                raise ValueError('Missing response_series.npz')
            rows[case], directories[case] = row, directory
        except (KeyError, ValueError) as exc:
            invalid.append({'case': case, 'error': str(exc)})
    pairs, curves = [], {}
    controls = [r for r in rows.values() if r['kind'] == 'control']
    for item in manifest['cases']:
        if item['kind'] != 'pulse' or item['case'] not in rows:
            continue
        case, control_case = item['case'], item['paired_control']
        pulse = rows[case]
        base = {key: item[key] for key in ('case', 'tm_over_tc', 'replicate')}
        base['initialization'] = pulse['initialization']
        base.update({'paired_control': control_case,
                     'pulse_duration_tc': pulse['pulse_duration_steps'] / TC})
        if control_case not in rows:
            base.update({'classification': 'missing_control', 'fits': {}, 'flags': []})
            pairs.append(base)
            continue
        alternate = next((r for r in controls
                          if r['tm_over_tc'] == pulse['tm_over_tc']
                          and r['replicate'] == pulse['replicate']
                          and r['initialization'] != pulse['initialization']), None)
        try:
            with np.load(directories[case] / 'response_series.npz') as ps, \
                    np.load(directories[control_case] / 'response_series.npz') as cs:
                result, response = pair_response(pulse, rows[control_case], ps, cs, alternate, allow_drift)
            base.update(result)
            if response is not None:
                curves[case] = response
        except (KeyError, ValueError) as exc:
            base.update({'classification': 'invalid_pair', 'fits': {}, 'flags': [str(exc)]})
        pairs.append(base)
    groups = []
    keys = sorted({(r['tm_over_tc'], str(r['initialization']), r['pulse_duration_tc']) for r in pairs})
    for tm, initialization, duration in keys:
        subset = [r for r in pairs if (r['tm_over_tc'], str(r['initialization']), r['pulse_duration_tc'])
                  == (tm, initialization, duration)]
        group = {'tm_over_tc': tm, 'initialization': initialization,
                 'pulse_duration_tc': duration, 'available_seed_count': len(subset),
                 'classification_counts': {label: sum(r['classification'] == label for r in subset)
                                           for label in sorted({r['classification'] for r in subset})},
                 'relaxation': {}}
        for metric in ('chi_mean', 'm_mean'):
            estimates = [{'replicate': r['replicate'], 'tau_tc': r['fits'][metric]['tau_tc']}
                         for r in subset if r['fits'].get(metric, {}).get('status') == 'resolved']
            values = [e['tau_tc'] for e in estimates]
            group['relaxation'][metric] = {'seed_estimates': estimates,
                                          'median_tau_tc': float(np.median(values)) if values else None,
                                          'range_tau_tc': [min(values), max(values)] if values else None}
        groups.append(group)
    consistency = []
    for tm, initialization in sorted({(g['tm_over_tc'], g['initialization']) for g in groups}):
        record = {'tm_over_tc': tm, 'initialization': initialization, 'metrics': {}}
        for metric in ('chi_mean', 'm_mean'):
            values = [{'pulse_duration_tc': g['pulse_duration_tc'],
                       'median_tau_tc': g['relaxation'][metric]['median_tau_tc'],
                       'resolved_seeds': len(g['relaxation'][metric]['seed_estimates'])}
                      for g in groups if g['tm_over_tc'] == tm and g['initialization'] == initialization
                      and g['relaxation'][metric]['median_tau_tc'] is not None]
            ratio = max(v['median_tau_tc'] for v in values) / min(v['median_tau_tc'] for v in values) if len(values) > 1 else None
            record['metrics'][metric] = {'by_duration': values, 'maximum_to_minimum_ratio': ratio,
                                         'interpretation': 'Descriptive duration sensitivity; no iid time-sample significance test.'}
        consistency.append(record)
    summary = {'expected_cases': len(manifest['cases']), 'available_cases': len(rows),
               'allow_drift': allow_drift,
               'missing': missing, 'invalid': invalid, 'cases': rows, 'pairs': pairs,
               'groups': groups, 'duration_consistency': consistency, 'criteria': CRITERIA,
               'interpretation': ('Fixed-age paired responses; baseline/velocity drift is descriptive. '
                    'Fit decay only when the paired late residual is consistent with zero. '
                    'Decay times are not proof of stationary-state stability. Seed ranges are not confidence intervals.'
                    if allow_drift else 'Finite-horizon basin and recovery evidence. No relaxation fit is assigned '
                    'to switched, drifting or unresolved responses. Seed ranges are not confidence intervals.')}
    out.mkdir(parents=True, exist_ok=True)
    for case, curve in curves.items():
        target = out / 'paired_series' / case
        target.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(target / 'paired_response.npz', **curve)
    write_json(out / 'pulse_summary.json', summary)
    if make_figures:
        figures(groups, pairs, curves, out)
    (out / 'README.md').write_text(
        '# Threshold-pulse response analysis\n\n'
        f'Available {len(rows)}/{len(manifest["cases"])} cases; '
        f'{len(missing)} missing and {len(invalid)} invalid.\n\n'
        'All pulse/control differences use matching absolute simulation steps. '
        f'Baseline screening uses the final {CRITERIA["baseline_windows_tc"]["far_from_transition"]:g} '
        f'tau_c before the pulse, or {CRITERIA["baseline_windows_tc"]["near_transition"]:g} tau_c '
        'near the transition. Frozen preparation is excluded. The entire feedback-on '
        'waiting window is also reported descriptively, but its initial transients '
        'do not gate baseline usability. '
        + ('Fixed-age mode treats all drift diagnostics as descriptive. Fits measure paired-response '
           'decay without requiring a stationary baseline; persistent late differences are kept without '
           'a forced zero-offset fit. These times alone do not establish stationary-state stability.\n\n'
           if allow_drift else 'Late-baseline and control drift make a response unsuitable for a relaxation fit. '
           'Terminal windows classify return, branch switching or incomplete recovery before fitting.\n\n') +
        'The late fit is A exp(-(t-t_start)/tau), with no constant offset. '
        'Initial continued excursions are excluded; all trial start times and fit quality '
        'are saved in pulse_summary.json. A fit must resolve decay above the noise floor '
        'over at least 1.5 fitted time constants and give stable tau for two late starts. '
        'Each seed is retained. Reported seed ranges and duration sensitivity are descriptive; '
        'time samples are correlated and four seeds do not establish precise switching probabilities.\n')
    print(json.dumps({'available': len(rows), 'expected': len(manifest['cases']),
                      'pairs': len(pairs), 'invalid': len(invalid)}))
    return summary


def figures(groups, pairs, curves, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'savefig.dpi': 180})
    durations = sorted({g['pulse_duration_tc'] for g in groups})
    palette = ['#63B8C6', '#8C6BB1', '#E8A0A6', '#68798C']
    colors = {d: palette[i % len(palette)] for i, d in enumerate(durations)}
    markers = {'patch': 's', '0': 'o', '1': '^'}
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7), constrained_layout=True)
    for ax, metric, label in zip(axes, ('chi_mean', 'm_mean'), (r'$\chi$', '$m$')):
        count = 0
        for g in groups:
            r = g['relaxation'][metric]
            if r['median_tau_tc'] is None:
                continue
            y = r['median_tau_tc']
            lo, hi = r['range_tau_tc']
            ax.errorbar(g['tm_over_tc'], y, yerr=[[y - lo], [hi - y]],
                        fmt=markers.get(str(g['initialization']), 'o'),
                        color=colors[g['pulse_duration_tc']], ms=5, lw=1., capsize=2)
            count += 1
        ax.set(xlabel=r'$\tau_m/\tau_c$', ylabel=r'$\tau_{\rm rel}/\tau_c$',
               title=label + ' relaxation')
        if count:
            ax.set_yscale('log')
        else:
            ax.text(.5, .5, 'No resolved exponential fits', transform=ax.transAxes,
                    ha='center', va='center', color='#68798C')
    handles = [Line2D([], [], color=colors[d], marker='o', linestyle='none',
                      label=rf'$t_p/\tau_c={d:.3g}$') for d in durations]
    handles += [Line2D([], [], color='#444444', marker=m, linestyle='none', label='init ' + a)
                for a, m in markers.items()]
    fig.legend(handles=handles, loc='outside lower center', ncol=4, frameon=False)
    for suffix in ('png', 'pdf'):
        fig.savefig(out / ('relaxation_times.' + suffix))
    plt.close(fig)
    available = [r for r in pairs if r['case'] in curves]
    if not available:
        return
    chosen = min(available, key=lambda r: (r['initialization'] != '1', abs(r['tm_over_tc'] - 8.)))
    selected = [r for r in available if r['tm_over_tc'] == chosen['tm_over_tc']
                and r['initialization'] == chosen['initialization']]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7), constrained_layout=True)
    for ax, metric, label in zip(axes, ('chi_mean', 'm_mean'), (r'$\Delta\langle\chi\rangle$', r'$\Delta\langle m\rangle$')):
        for duration in durations:
            rr = [r for r in selected if r['pulse_duration_tc'] == duration]
            binned = []
            for row in rr:
                curve = curves[row['case']]
                t, y = bin_series(curve['time_after_pulse_tc'], curve['delta_' + metric], .25)
                mask = t <= 50.
                ax.plot(t[mask], y[mask], color=colors[duration], alpha=.25, lw=.6)
                binned.append((t, y))
            if binned:
                # Group members share duration and clock, so bin arrays align.
                length = min(len(t) for t, _ in binned)
                t = binned[0][0][:length]
                y = np.mean([values[:length] for _, values in binned], axis=0)
                mask = t <= 50.
                ax.plot(t[mask], y[mask], color=colors[duration], lw=1.5,
                        label=rf'$t_p/\tau_c={duration:.3g}$')
        ax.axhline(0, color='#A0A0A0', lw=.6)
        ax.set(xlabel=r'$(t-t_{\rm off})/\tau_c$', ylabel=label)
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle(f'Representative paired responses: tau_m/tau_c={chosen["tm_over_tc"]:g}, '
                 f'init {chosen["initialization"]}; thin lines: seeds')
    for suffix in ('png', 'pdf'):
        fig.savefig(out / ('representative_responses.' + suffix))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--summary', action='store_true')
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--no-figures', action='store_true')
    parser.add_argument('--allow-drift', action='store_true',
                        help='Analyze fixed-age paired responses without stationarity/velocity vetoes')
    args = parser.parse_args()
    if args.summary:
        if args.manifest is None:
            parser.error('--summary requires --manifest')
        aggregate(args.input, args.output, args.manifest, not args.no_figures, args.allow_drift)
    else:
        reduce_case(args.input, args.output)


if __name__ == '__main__':
    main()
