#!/usr/bin/env python3
"""Audit and summarize the three-seed, 10%-p_c poster pulse campaign.

Per case: RAW_CASE OUT_CASE --manifest MANIFEST
Summary:  RESULTS_ROOT OUT_SUMMARY --manifest MANIFEST --summary

The per-case reducer and exponential-fit criteria are inherited unchanged from
cw_pmem_pulse_analysis.py. The primary estimate fits the signed response averaged
over ALL THREE seeds. Drift flags are retained, never used to choose seeds.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

import cw_pmem_pulse_analysis as legacy


METRICS = ('chi_mean', 'm_mean')
ALIASES = {'model': 'model_name', 'bc': 'BC', 'initial-order': 'init_order'}
# These disabled tracer controls are intentionally absent from serialize_params.
UNSERIALIZED_DISABLED = {'ntracer', 'tracer-count'}
INTEGER_PARAMETERS = {'LX', 'LY', 'BC', 'nsteps', 'nstart', 'ninfo', 'nsubsteps',
                      'seed', 'chi_seed', 'chi_freeze_steps', 'mem_freeze_steps',
                      'pmem_pulse_start', 'pmem_pulse_steps', 'nresponse', 'nvideo',
                      'video_stride', 'nvideo_dense', 'video_start', 'video_dense_start',
                      'video_dense_end', 'open_loop', 'frame_light', 'switch_sign'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def equal(actual, expected):
    """Runtime JSON uses six significant digits for floating-point values."""
    try:
        return bool(np.isclose(float(actual), float(expected), rtol=6e-6, atol=1e-12))
    except (TypeError, ValueError):
        return str(actual) == str(expected)


def load_manifest(path):
    manifest = read_json(path)
    fixed = manifest['fixed']
    if manifest['replicates'] != 3:
        raise ValueError('This campaign requires exactly three independent seeds per group')
    required = {'L': 256, 'mc': .2287, 'pmem': .016838,
                'pmem_pulse_value': .0151542, 'tau_chi': 202.3,
                'r': .3, 'tau_c': legacy.TC}
    for key, value in required.items():
        if key not in fixed or not equal(fixed[key], value):
            raise ValueError(f'Unexpected campaign fixed parameter {key}')
    if not equal(fixed['pmem_pulse_value'], .9 * fixed['pmem']):
        raise ValueError('The pulse must lower p_c by 10%')
    if manifest['tp_grid'] != [3] or not manifest['tm_grid']:
        raise ValueError('Expected a nonempty memory grid and t_p/tau_c=3')
    names = [item['case'] for item in manifest['cases']]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate cases in manifest')
    index = {item['case']: item for item in manifest['cases']}
    expected = {(float(g), init, rep) for g in manifest['tm_grid']
                for init in ('0', '1') for rep in (1, 2, 3)}
    for kind in ('control', 'pulse'):
        members = [item for item in index.values() if item['kind'] == kind]
        actual = [(float(item['tm_over_tc']), str(item['initialization']),
                   int(item['replicate'])) for item in members]
        if len(actual) != len(expected) or set(actual) != expected:
            raise ValueError(f'Incomplete or duplicated manifest grid for {kind}')
    for item in index.values():
        name = Path(item['case'])
        if name.is_absolute() or '..' in name.parts:
            raise ValueError('Case must be a relative path below the campaign')
        if item['kind'] == 'pulse':
            control = index[item['paired_control']]
            if control['kind'] != 'control':
                raise ValueError('paired_control is not a control')
            for key in ('tm_over_tc', 'initialization', 'replicate', 'seed', 'chi_seed',
                        'preparation_steps', 'feedback_wait_steps', 'pulse_start_steps', 'nsteps'):
                if str(item[key]) != str(control[key]):
                    raise ValueError(f'Manifest pair mismatch: {key}')
            if abs(item['pulse_duration_steps'] / legacy.TC - 3.) > 1 / legacy.TC:
                raise ValueError('Wrong pulse duration')
        elif item['pulse_duration_steps'] != 0:
            raise ValueError('Control has a nonzero pulse duration')
        if abs(item['feedback_wait_steps'] / legacy.TC - 2000.) > 1 / legacy.TC:
            raise ValueError('All groups must share a 2000 tau_c feedback wait')
        if item['pulse_start_steps'] != item['preparation_steps'] + item['feedback_wait_steps']:
            raise ValueError('Inconsistent feedback-release/pulse-start clock')
        if item['nsteps'] - item['pulse_start_steps'] - item['pulse_duration_steps'] < 999.99 * legacy.TC:
            raise ValueError('Less than 1000 tau_c recovery (within step rounding)')
    for g in manifest['tm_grid']:
        for init in ('0', '1'):
            members = [item for item in index.values() if item['kind'] == 'pulse'
                       and item['tm_over_tc'] == g and str(item['initialization']) == init]
            if len({item['seed'] for item in members}) != 3:
                raise ValueError('Seeds within a group must be distinct')
    return manifest


def audit_parameters(parameters, item, fixed):
    expected = item['expected_parameters']
    if not expected:
        raise ValueError('Missing expected runtime parameters')
    for key, value in expected.items():
        actual_key = ALIASES.get(key, key.replace('-', '_'))
        if key in UNSERIALIZED_DISABLED and actual_key not in parameters and equal(value, 0):
            continue
        if actual_key not in parameters:
            raise ValueError(f'Missing serialized parameter: {actual_key}')
        matches = (float(parameters[actual_key]) == float(value) if actual_key in INTEGER_PARAMETERS
                   else equal(parameters[actual_key], value))
        if not matches:
            raise ValueError(f'Runtime mismatch {actual_key}: {parameters[actual_key]} != {value}')
    checks = {'LX': fixed['L'], 'LY': fixed['L'], 'mc': fixed['mc'],
              'pmem': fixed['pmem'], 'pmem_pulse_value': .9 * fixed['pmem'],
              'tau_chi': fixed['tau_chi'], 'zeta0_frac': fixed['r'],
              'tau_m': item['tm_over_tc'] * legacy.TC, 'chi_config': 'uniform',
              'chi0': item['initialization'], 'seed': item['seed'],
              'chi_seed': item['chi_seed'], 'nsteps': item['nsteps'],
              'chi_freeze_steps': item['preparation_steps'],
              'pmem_pulse_start': item['pulse_start_steps'],
              'pmem_pulse_steps': item['pulse_duration_steps'], 'open_loop': 0}
    for key, value in checks.items():
        matches = key in parameters and (float(parameters[key]) == float(value)
                  if key in INTEGER_PARAMETERS else equal(parameters[key], value))
        if not matches:
            raise ValueError(f'Runtime identity/protocol mismatch: {key}')


def reduce_case(root, out, manifest_path):
    manifest = load_manifest(manifest_path)
    matches = [item for item in manifest['cases']
               if root.as_posix().rstrip('/').endswith('/' + item['case'])]
    if len(matches) != 1:
        raise ValueError(f'Cannot identify exactly one manifest case for {root}')
    item = matches[0]
    card = manifest_path.parent / item['case'] / 'run.dat'
    if sha(card) != item['run_dat_sha256']:
        raise ValueError('Tracked run.dat hash differs from manifest')
    parameters = legacy.parameters(root)
    audit_parameters(parameters, item, manifest['fixed'])
    legacy.reduce_case(root, out)
    provenance = {
        'case': item['case'], 'manifest_sha256': sha(manifest_path),
        'run_dat_sha256': sha(card), 'parameters_sha256': sha(root / 'parameters.json'),
        'response_csv_sha256': sha(root / 'response.csv'),
        'pulse_case_sha256': sha(out / 'pulse_case.json'),
        'response_series_sha256': sha(out / 'response_series.npz'),
        'reducer_sha256': sha(Path(legacy.__file__)), 'analysis_sha256': sha(Path(__file__)),
        'parameters': parameters,
        'unserialized_disabled_controls': {key: item['expected_parameters'][key]
                                           for key in UNSERIALIZED_DISABLED
                                           if key in item['expected_parameters']},
        'relative_threshold_decrease': 1 - parameters['pmem_pulse_value'] / parameters['pmem'],
    }
    write_json(out / 'poster_pulse_provenance.json', provenance)


def locate(root, case):
    for directory in (root / case, root / case / 'analysis'):
        if (directory / 'pulse_case.json').is_file():
            return directory
    return None


def load_case(directory, item, manifest, manifest_path):
    row = read_json(directory / 'pulse_case.json')
    provenance = read_json(directory / 'poster_pulse_provenance.json')
    expected_hashes = {'manifest_sha256': sha(manifest_path),
                       'run_dat_sha256': item['run_dat_sha256'],
                       'pulse_case_sha256': sha(directory / 'pulse_case.json'),
                       'response_series_sha256': sha(directory / 'response_series.npz'),
                       'reducer_sha256': sha(Path(legacy.__file__)),
                       'analysis_sha256': sha(Path(__file__))}
    for key, expected in expected_hashes.items():
        if provenance.get(key) != expected:
            raise ValueError(f'Invalid provenance: {key}')
    if provenance['case'] != item['case'] or row['case_name'] != Path(item['case']).name:
        raise ValueError('Mismatched case identity')
    if not row['complete']:
        raise ValueError('Incomplete case')
    audit_parameters(provenance['parameters'], item, manifest['fixed'])
    for key in ('kind', 'initialization', 'replicate', 'seed', 'nsteps',
                'preparation_steps', 'pulse_start_steps', 'pulse_duration_steps'):
        if str(row[key]) != str(item[key]):
            raise ValueError(f'Reduced metadata mismatch: {key}')
    for key in ('L', 'mc', 'pmem', 'pmem_pulse_value', 'tau_chi'):
        if not equal(row[key], manifest['fixed'][key]):
            raise ValueError(f'Reduced fixed parameter mismatch: {key}')
    if not equal(row['tm_over_tc'], item['tm_over_tc']) or not equal(row['tau_c_reference'], legacy.TC):
        raise ValueError('Reduced clock mismatch')
    return row, provenance


def analyze_mean(time, seed_response, control_blocks, metric):
    """No offset subtraction, normalization, or single-seed fit-status selection."""
    values = np.asarray(seed_response, dtype=float)
    blocks = np.asarray(control_blocks, dtype=float)
    if values.shape != (3, len(time)) or blocks.shape != (3, 5):
        raise ValueError('The ensemble requires all three synchronized seed responses')
    if not np.isfinite(values).all() or not np.isfinite(blocks).all():
        raise ValueError('Non-finite ensemble member')
    mean = np.mean(values, axis=0)
    sem = np.std(values, axis=0, ddof=1) / np.sqrt(3.)
    tail = time >= time[-1] - min(legacy.CRITERIA['terminal_window_max_tc'], time[-1] / 2.)
    residual = float(np.mean(mean[tail]))
    noise = float(np.std(np.mean(blocks, axis=0)))
    tolerance = max(legacy.CRITERIA['return_absolute_tolerances'][metric],
                    legacy.CRITERIA['return_noise_multiplier'] * noise)
    returned = abs(residual) <= tolerance
    fit = legacy.exponential_fit(time, mean) if returned else {
        'status': 'persistent_mean_residual', 'tau_tc': None, 'candidates': []}
    terminal_seeds = np.mean(values[:, tail], axis=1)
    result = {'terminal_mean': residual, 'return_tolerance': tolerance,
              'control_mean_block_std': noise, 'mean_returned_within_tolerance': bool(returned),
              'terminal_seed_values': terminal_seeds.tolist(),
              'terminal_seed_sem': float(np.std(terminal_seeds, ddof=1) / np.sqrt(3.)),
              'fit': fit}
    return result, mean, sem


def validate_pair_clock(curve, pulse):
    time, step = curve['time_after_pulse_tc'], curve['step']
    if len(time) < 30 or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError('Invalid paired recovery clock')
    if not np.array_equal(time, (step - pulse['pulse_end_steps']) / legacy.TC):
        raise ValueError('Incorrect time origin or time units')
    if step[-1] != pulse['nsteps'] or step[0] < pulse['pulse_end_steps']:
        raise ValueError('Incomplete recovery coverage')
    for metric in METRICS:
        values = curve['delta_' + metric]
        if values.shape != time.shape or not np.isfinite(values).all() or np.max(np.abs(values)) > 1 + 1e-8:
            raise ValueError('Invalid signed paired response: ' + metric)


def group_name(tm, init, tp):
    return f'tm{tm:g}_init{init}_tp{tp:g}'.replace('.', 'p')


def summarize(root, out, manifest_path, make_figures=True):
    manifest = load_manifest(manifest_path)
    out.mkdir(parents=True, exist_ok=True)
    rows, provenance, directories, coverage = {}, {}, {}, []
    for item in manifest['cases']:
        case = item['case']
        record = {'case': case, 'kind': item['kind'], 'tm_over_tc': item['tm_over_tc'],
                  'initialization': str(item['initialization']), 'replicate': item['replicate'],
                  'seed': item['seed'], 'status': 'missing', 'error': None}
        directory = locate(root, case)
        if directory is not None:
            try:
                row, proof = load_case(directory, item, manifest, manifest_path)
                rows[case], provenance[case], directories[case] = row, proof, directory
                record['status'] = 'complete'
                record['baseline_flags'] = row['baseline']['flags']
                record['control_post_flags'] = row['control_post_flags']
                record['terminal_flags'] = row['terminal']['flags']
            except (KeyError, ValueError, OSError) as exc:
                record.update(status='invalid', error=str(exc))
        coverage.append(record)
    write_json(out / 'coverage.json', {'expected': len(coverage), 'cases': coverage})
    with (out / 'coverage.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, ['case', 'kind', 'tm_over_tc', 'initialization',
                                        'replicate', 'seed', 'status', 'error'], extrasaction='ignore')
        writer.writeheader()
        writer.writerows(coverage)

    expected_groups = defaultdict(list)
    for item in manifest['cases']:
        if item['kind'] == 'pulse':
            expected_groups[(float(item['tm_over_tc']), str(item['initialization']),
                             float(item['tp_over_tc']))].append(item)
    groups, flat, pair_reports = [], [], []
    curve_root = out / 'curves'
    curve_root.mkdir(exist_ok=True)
    for (tm, init, tp), members in sorted(expected_groups.items()):
        members.sort(key=lambda item: item['replicate'])
        name = group_name(tm, init, tp)
        group = {'id': name, 'tm_over_tc': tm, 'initialization': init, 'tp_over_tc': tp,
                 'seeds': [item['seed'] for item in members], 'expected_seeds': 3,
                 'available_seeds': 0, 'status': 'incomplete_group', 'errors': [],
                 'members': [], 'metrics': {}}
        all_values = {metric: [] for metric in METRICS}
        control_blocks = {metric: [] for metric in METRICS}
        reference_step, reference_time = None, None
        for item in members:
            case, control_case = item['case'], item['paired_control']
            if case not in rows or control_case not in rows:
                group['errors'].append({'case': case, 'error': 'missing_or_invalid_pair_member'})
                continue
            try:
                pulse, control = rows[case], rows[control_case]
                pp, cp = provenance[case]['parameters'], provenance[control_case]['parameters']
                for key in pp:
                    # The scheduled pulse duration is the only runtime difference.
                    if key != 'pmem_pulse_steps' and pp[key] != cp.get(key):
                        raise ValueError(f'Runtime pulse/control mismatch: {key}')
                with np.load(directories[case] / 'response_series.npz') as ps, \
                        np.load(directories[control_case] / 'response_series.npz') as cs:
                    report, curve = legacy.pair_response(pulse, control, ps, cs, allow_drift=True)
                if curve is None or report['flags']:
                    raise ValueError('Invalid pair: ' + ', '.join(report['flags']))
                validate_pair_clock(curve, pulse)
                if reference_step is None:
                    reference_step = curve['step'].copy()
                    reference_time = curve['time_after_pulse_tc'].copy()
                elif not (np.array_equal(curve['step'], reference_step)
                          and np.array_equal(curve['time_after_pulse_tc'], reference_time)):
                    raise ValueError('Seed curves do not have identical absolute steps and time grids')
                target = out / 'paired_series' / case
                target.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(target / 'paired_response.npz', **curve)
                pair_reports.append({'case': case, 'paired_control': control_case, **report})
                group['members'].append({'case': case, 'paired_control': control_case,
                    'seed': item['seed'], 'replicate': item['replicate'],
                    'classification': report['classification'],
                    'diagnostic_flags': report['diagnostic_flags'],
                    'baseline_stationary_by_legacy_checks': report['baseline_stationary_by_legacy_checks'],
                    'pulse_baseline_flags': pulse['baseline']['flags'],
                    'pulse_terminal_flags': pulse['terminal']['flags'],
                    'control_baseline_flags': control['baseline']['flags'],
                    'control_post_flags': control['control_post_flags'],
                    'control_terminal_flags': control['terminal']['flags']})
                for metric in METRICS:
                    all_values[metric].append(curve['delta_' + metric])
                    control_blocks[metric].append(control['terminal']['metrics'][metric]['blocks'])
            except (KeyError, ValueError, OSError) as exc:
                group['errors'].append({'case': case, 'error': str(exc)})
        group['available_seeds'] = len(group['members'])
        if len(group['members']) == 3 and not group['errors']:
            group.update(status='complete', n_seeds=3,
                         feedback_wait_tc=members[0]['feedback_wait_steps'] / legacy.TC,
                         recovery_horizon_tc=float(reference_time[-1]))
            exported = {'step': reference_step, 'time_after_pulse_tc': reference_time,
                        'seed': np.asarray(group['seeds'])}
            for metric in METRICS:
                values = np.asarray(all_values[metric])
                result, mean, sem = analyze_mean(reference_time, values, control_blocks[metric], metric)
                group['metrics'][metric] = result
                exported[metric + '_seed_response'] = values
                exported[metric + '_mean'] = mean
                exported[metric + '_sem'] = sem
                exported[metric + '_control_terminal_blocks'] = np.asarray(control_blocks[metric])
            np.savez_compressed(curve_root / (name + '.npz'), **exported)
            group['curve_file'] = 'curves/' + name + '.npz'
            group['curve_sha256'] = sha(out / group['curve_file'])
        else:
            for metric in METRICS:
                group['metrics'][metric] = {'fit': {'status': 'incomplete_group', 'tau_tc': None}}
        for metric in METRICS:
            result = group['metrics'][metric]
            flat.append({'id': name, 'tm_over_tc': tm, 'initialization': init, 'tp_over_tc': tp,
                         'metric': metric, 'n_seeds': group['available_seeds'],
                         'status': result['fit']['status'], 'tau_tc': result['fit'].get('tau_tc'),
                         'terminal_mean': result.get('terminal_mean'),
                         'terminal_sem': result.get('terminal_seed_sem'),
                         'return_tolerance': result.get('return_tolerance'),
                         'drift_flagged_members': sum(bool(member['diagnostic_flags']) for member in group['members'])})
        groups.append(group)
        print(json.dumps({'group': name, 'available_seeds': group['available_seeds'],
                          'fits': {metric: group['metrics'][metric]['fit']['status'] for metric in METRICS}}), flush=True)
    result = {'manifest': str(manifest_path), 'manifest_sha256': sha(manifest_path),
              'analysis_sha256': sha(Path(__file__)), 'fit_script_sha256': sha(Path(legacy.__file__)),
              'fixed': manifest['fixed'], 'expected_cases': len(manifest['cases']),
              'available_cases': len(rows), 'expected_groups': len(expected_groups),
              'complete_groups': sum(group['status'] == 'complete' for group in groups),
              'independent_seeds_per_group': 3, 'criteria': legacy.CRITERIA,
              'definition': 'R_X(t) = mean over all 3 seeds of [X_pulse(step)-X_control(step)], t=(step-pulse_end)/tau_c.',
              'averaging': 'Signed responses averaged before fitting; no fit-based member selection, interpolation, normalization, or offset removal.',
              'terminal_check': 'Independent per observable; |terminal mean| <= max(legacy absolute tolerance, 2 std of five synchronized ensemble-control terminal block means).',
              'fit_model': 'A exp[-(t-t_start)/tau]; fixed zero offset; inherited fit windows and acceptance criteria.',
              'uncertainty': 'Pointwise sample SD/sqrt(3) across independent seeds; SEM is not a confidence band. No relaxation-time confidence interval or leave-one-out interval assigned.',
              'interpretation': 'Fixed-age paired response; all baseline/control/pulse drift flags retained without seed selection. A fitted decay does not establish stationary-state stability or linear response.',
              'coverage': coverage, 'provenance': provenance, 'pairs': pair_reports,
              'groups': groups, 'fit_counts': {metric: dict(Counter(group['metrics'][metric]['fit']['status'] for group in groups)) for metric in METRICS}}
    write_json(out / 'ensemble_results.json', result)
    with (out / 'ensemble_fits.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    write_readme(out, result)
    if make_figures:
        figures(out, result)
    complete = result['complete_groups'] == result['expected_groups']
    print(json.dumps({'complete': complete, 'groups': result['complete_groups'],
                      'expected_groups': result['expected_groups'], 'fit_counts': result['fit_counts']}), flush=True)
    return complete


def write_readme(out, result):
    (out / 'README.md').write_text(
        '# Dense 10% threshold-pulse response\n\n'
        f"Coverage: {result['available_cases']}/{result['expected_cases']} cases; "
        f"{result['complete_groups']}/{result['expected_groups']} complete three-seed groups.\n\n"
        'L256, mc=0.2287, pc=0.016838, activity floor r=0.3, tau_chi=202.3 steps. '
        'pc is temporarily lowered to 0.9 pc, for 3 tau_c (integer-step rounding), '
        'after 2000 tau_c feedback-on waiting. Both uniform initializations are shown '
        'separately, with 1000 tau_c recovery. tau_c=674.3290333006435 steps.\n\n'
        'At each memory time and initialization, the signed pulse-minus-control curves '
        'are averaged over all three seeds BEFORE fitting. Controls have the same seed '
        'and deterministic prehistory. Pairing uses identical absolute simulation steps, '
        'without interpolation; all three paired grids must match exactly. No member is '
        'selected by fit quality, stationarity, sign or amplitude. Missing or invalid '
        'members leave the whole group unestimated (never zero or infinity).\n\n'
        'The operational terminal return test is independent for chi and m: the late mean '
        'must be within max(original absolute tolerance, twice the standard deviation of '
        'five synchronized ensemble-control terminal block means). This is not a significance '
        'test. Persistent residuals are recorded without subtracting an offset. Fits use '
        'the unchanged cw_pmem_pulse_analysis.exponential_fit and 0.1 tau_c bins. '
        'All fit-window candidates and rejection reasons are retained.\n\n'
        'Raw arrays, individual seed curves, their mean and pointwise SEM are saved in '
        'curves/*.npz. SEM is sample standard deviation / sqrt(3), not a 95% confidence '
        'band. Figures use 1 tau_c bins only for display; no extra smoothing or centering. '
        'The curve atlas includes every available complete group and both observables, '
        'with individual seeds and the mean plus/minus SEM. Relaxation scatter contains '
        'only resolved fits; coverage and all failures remain in the CSV/JSON and status figure.\n\n'
        'These are fixed-age paired responses. Baseline, control and pulse-terminal drift '
        'flags remain in the report and do not exclude seeds. A three-seed mean can hide '
        'individual persistent responses; inspect individual curves. A finite 10% pulse '
        'does not by itself verify linear response, a stationary-state restoring rate, '
        'a critical exponent or divergence. No time-constant uncertainty is inferred '
        'from correlated time samples.\n\n'
        'Provenance includes manifest and run.dat hashes, all serialized runtime parameters, '
        'parameters.json and response.csv hashes, reduced-data hashes, and reducer/analysis hashes. '
        'Runtime p_c and pulse p_c are audited against the manifest, as are all serialized '
        'expected parameters. The summary exits unsuccessfully if a group is incomplete, '
        'while still writing coverage and estimates for complete groups.\n')


def figures(out, report):
    """Analysis summaries; fit quantities and display binning remain separate."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': False, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    colors, markers = {'0': '#2471A3', '1': '#B03A2E'}, {'0': 'o', '1': 's'}
    groups = report['groups']

    def save(fig, name):
        for extension in ('pdf', 'png'):
            fig.savefig(out / (name + '.' + extension), dpi=200, facecolor='white')
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.1))
    fig.subplots_adjust(left=.085, right=.98, bottom=.25, top=.91, wspace=.28)
    for ax, metric, symbol in zip(axes, METRICS, (r'\chi', 'm')):
        count = 0
        for group in groups:
            fit = group['metrics'][metric]['fit']
            if fit['status'] == 'resolved':
                init = group['initialization']
                ax.plot(group['tm_over_tc'], fit['tau_tc'], marker=markers[init],
                        color=colors[init], markerfacecolor='white', linestyle='none', ms=5)
                count += 1
        ax.set(xlabel=r'$\tau_m/\tau_c$', ylabel=r'$\tau_{\mathrm{rel},' + symbol + r'}/\tau_c$')
        if count:
            ax.set_yscale('log')
        else:
            ax.text(.5, .5, 'No resolved mean-response fits', transform=ax.transAxes, ha='center')
    handles = [Line2D([], [], color=colors[init], marker=markers[init],
                      mfc='white', linestyle='none', label='initialization ' + init) for init in ('0', '1')]
    fig.legend(handles=handles, loc='upper center', ncol=2, frameon=False)
    fig.text(.085, .095, 'Each point fits the signed mean of all 3 seeds; no offset or fitted-time averaging.', fontsize=9)
    fig.text(.085, .045, 'Missing points have no assigned relaxation time; see fit_status and coverage tables.', fontsize=9)
    save(fig, 'ensemble_relaxation')

    statuses = ['resolved', 'persistent_mean_residual', 'signal_too_weak',
                'decay_not_resolved', 'incomplete_group', 'other_unresolved_fit']
    status_colors = ['#2471A3', '#B03A2E', '#B2BABB', '#D68910', '#222222', '#7D3C98']
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.subplots_adjust(left=.09, right=.98, bottom=.35, top=.85, wspace=.28)
    for ax, metric, title in zip(axes, METRICS, ('chi', 'memory')):
        for group in groups:
            status = group['metrics'][metric]['fit']['status']
            i = statuses.index(status) if status in statuses else len(statuses) - 1
            ax.scatter(group['tm_over_tc'], int(group['initialization']), c=status_colors[i], s=35,
                       marker=markers[group['initialization']])
        ax.set(xlabel=r'$\tau_m/\tau_c$', yticks=[0, 1], ylim=(-.4, 1.4), title=title)
        ax.set_ylabel('initialization')
    fig.legend(handles=[Line2D([], [], color=c, marker='o', linestyle='none', label=s.replace('_', ' '))
                        for s, c in zip(statuses, status_colors)], loc='lower center',
               ncol=2, frameon=False, fontsize=9)
    save(fig, 'fit_status')

    complete = [group for group in groups if group['status'] == 'complete']
    if not complete:
        return
    with PdfPages(out / 'ensemble_response_atlas.pdf') as pdf:
        for metric, symbol in zip(METRICS, (r'\chi', 'm')):
            for start in range(0, len(complete), 8):
                page = complete[start:start + 8]
                nrows = (len(page) + 1) // 2
                height = 2.7 * nrows + 1.2
                fig, axes = plt.subplots(nrows, 2, figsize=(11.7, height), squeeze=False)
                fig.subplots_adjust(left=.085, right=.98, bottom=1.0 / height,
                                    top=1 - .8 / height, hspace=.62, wspace=.27)
                for ax, group in zip(axes.flat, page):
                    with np.load(out / group['curve_file']) as data:
                        time = data['time_after_pulse_tc']
                        values = data[metric + '_seed_response']
                        bt, _ = legacy.bin_series(time, values[0], 1.)
                        bv = np.asarray([legacy.bin_series(time, seed_values, 1.)[1]
                                         for seed_values in values])
                    mean, sem = bv.mean(axis=0), bv.std(axis=0, ddof=1) / np.sqrt(3.)
                    color = colors[group['initialization']]
                    for values in bv:
                        ax.plot(bt, values, color=color, lw=.5, alpha=.25)
                    ax.fill_between(bt, mean - sem, mean + sem, color=color, alpha=.18, linewidth=0)
                    ax.plot(bt, mean, color=color, lw=1.4)
                    ax.axhline(0, color='#777777', lw=.6, ls=':')
                    fit = group['metrics'][metric]['fit']
                    if fit['status'] == 'resolved':
                        ft = np.linspace(fit['fit_start_tc'], fit['fit_end_tc'], 150)
                        fy = fit['amplitude_at_start'] * np.exp(-(ft - fit['fit_start_tc']) / fit['tau_tc'])
                        ax.plot(ft, fy, 'k--', lw=1)
                        label = f'tau_rel = {fit["tau_tc"]:.3g} tau_c'
                    else:
                        label = fit['status'].replace('_', ' ')
                    flagged = sum(bool(member['diagnostic_flags']) for member in group['members'])
                    ax.set(title=f'g = {group["tm_over_tc"]:g}; init {group["initialization"]}; drift flags {flagged}/3',
                           xlabel=r'$(t-t_{\mathrm{off}})/\tau_c$', ylabel=r'$R_' + symbol + '$',
                           xlim=(0, group['recovery_horizon_tc']))
                    ax.text(.98, .92, label, transform=ax.transAxes, ha='right', va='top', fontsize=8,
                            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .85})
                for ax in list(axes.flat)[len(page):]:
                    ax.set_visible(False)
                fig.suptitle(r'$R_' + symbol + r'$: signed mean over all 3 seeds; $\Delta p_c/p_c=-0.1$, $t_p=3\tau_c$',
                             y=1 - .14 / height)
                fig.text(.085, .36 / height, 'Thin: individual seeds. Solid/shading: mean +/- SEM. Dashed: accepted zero-offset fit.', fontsize=10)
                fig.text(.085, .16 / height, 'Display bins: 1 tau_c; fit bins: 0.1 tau_c. Drift flags do not exclude seeds. SEM is not a confidence interval.', fontsize=9)
                pdf.savefig(fig)
                if start == 0:
                    fig.savefig(out / f'response_atlas_{metric}_preview.png', dpi=160, facecolor='white')
                plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--summary', action='store_true')
    parser.add_argument('--no-figures', action='store_true')
    args = parser.parse_args()
    manifest = args.manifest.resolve()
    if args.summary:
        if not summarize(args.root, args.out, manifest, not args.no_figures):
            raise SystemExit(2)
    else:
        reduce_case(args.root.resolve(), args.out, manifest)


if __name__ == '__main__':
    main()
