#!/usr/bin/env python3
"""Generate matched open-loop distributions and dense 10% sensing pulses.

No simulator equations are changed. The source baseline is the mc=.2287 core
figure; fields that differ are experiment protocols or output controls.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import subprocess

REPO = Path(__file__).resolve().parents[1]
BASELINE = Path('cases/20260911/cw_mc02287_four_init/L256/mc0p2287_tm0p3_chi0_rep1/run.dat')
TC = 674.3290333006435
MC = .2287
PC = .016838
ACTIVITIES = (.3, .5, .7, 1.)
DISTRIBUTION_TM = (.3, 1., 3., 6., 8., 10., 20., 30.)
PULSE_TM = (1., 2., 3., 4., 5., 6., 6.5, 7., 7.25, 7.5, 7.7, 7.9,
            8.1, 8.3, 8.5, 8.75, 9., 9.5, 10., 11., 12., 14., 16., 18., 20.)
DIST_CAMPAIGN = '20260917/cw_poster_distributions'
PULSE_CAMPAIGN = '20260917/cw_poster_pulse10'
PHYSICS_KEYS = ('model', 'LX', 'LY', 'bc', 'Gamma', 'xi', 'CC', 'LL', 'tau',
                'rho', 'friction', 'zeta', 'zeta0-frac', 'angle', 'noise',
                'initial-order', 'director-config', 'defect-sep', 'Dbio',
                'tau-chi', 'chi-width', 'switch-sign', 'mc', 'pmem', 'pmem-width')


def parse_card(text):
    result = {}
    for line in text.splitlines():
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        key, value = line.split('=', 1)
        if key.strip() in result:
            raise ValueError(f'Duplicate parameter {key}')
        result[key.strip()] = value.strip()
    return result


def token(value):
    return f'{value:g}'.replace('.', 'p')


def source_provenance():
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=REPO, text=True).splitlines()
    hashes = {name: hashlib.sha256(subprocess.check_output(
        ['git', 'show', f'HEAD:{name}'], cwd=REPO)).hexdigest() for name in names}
    return dict(source_base_commit=commit, simulator_source_sha256=hashes,
                baseline_runcard=str(BASELINE),
                baseline_runcard_sha256=hashlib.sha256((REPO/BASELINE).read_bytes()).hexdigest())


def immutable_write(path, text):
    if path.exists() and path.read_text() != text:
        raise FileExistsError(f'Refusing to replace differing generated file: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_case(root, name, parameters, comments, base):
    values = {key: f'{value:.16g}' if isinstance(value, float) else str(value)
              for key, value in parameters.items()}
    if any(values[key] != base[key] for key in PHYSICS_KEYS):
        raise ValueError(f'Core physics differs: {name}')
    text = ''.join(f'# {line}\n' for line in comments) + ''.join(
        f'{key:24s} = {value}\n' for key, value in values.items())
    immutable_write(root/name/'run.dat', text)
    return dict(case=name, expected_parameters=values,
                run_dat_sha256=hashlib.sha256(text.encode()).hexdigest())


def common_manifest(campaign, base, provenance):
    return dict(campaign=campaign, **provenance, tau_c=TC, replicates=3,
                fixed=dict(L=256, mc=MC, pmem=PC, tau_chi=202.3, r=.3, tau_c=TC),
                reference_physics={k: base[k] for k in PHYSICS_KEYS},
                simulator_note='Use committed simulator only. Unrelated local model edits excluded.',
                seed_note='Three independent seeds per condition; seeds paired across parameter points. '
                          'Parameter points and spatial sites are not independent seed replicates.')


def distribution_campaign(parent, base, provenance):
    root = parent/DIST_CAMPAIGN
    rows = []
    ninfo = round(5*TC)
    sample_steps = math.ceil(600*TC/ninfo)*ninfo
    for tm in DISTRIBUTION_TM:
        steady_start = math.ceil(max(200., 20*tm)*TC/ninfo)*ninfo
        nsteps = steady_start+sample_steps
        for activity in ACTIVITIES:
            for rep in (1, 2, 3):
                p = dict(base)
                p.update({'open-loop': 1, 'zeta-open': activity*float(base['zeta']),
                          'seed': 191700+rep, 'chi-seed': 201700+rep,
                          'chi-config': 'uniform', 'chi0': .5, 'm0': .25,
                          'chi-freeze-steps': 0, 'mem-freeze-steps': 0,
                          'tau-m': tm*TC, 'nsteps': nsteps, 'ninfo': ninfo, 'nstart': 0,
                          'nvideo': 0, 'ntracer': 0, 'tracer-count': 0,
                          'nresponse': 0, 'frame-light': 1,
                          'pmem-pulse-steps': 0})
                name = f'L256/tm{token(tm)}_a{token(activity)}_rep{rep}'
                row = write_case(root, name, p, [DIST_CAMPAIGN,
                    'Uniform prescribed activity; memory and phenotype evolve without stress feedback.',
                    f'tau_c={TC:.16g}; tau_m/tau_c={tm:g}; a={activity:g}; pc={PC:g}.',
                    f'Steady distribution window [{steady_start},{nsteps}], cadence {ninfo} steps.',
                    'Native m and pressure share the same frames. Match mean and second moment; no MLE.'], base)
                row.update(activity=activity, tm_over_tc=tm, seed=191700+rep,
                           replicate=rep, steady_start=steady_start, nsteps=nsteps,
                           ninfo=ninfo, steady_frames=sample_steps//ninfo+1)
                rows.append(row)
    manifest = common_manifest(DIST_CAMPAIGN, base, provenance)
    manifest.update(activity_grid=list(ACTIVITIES), tm_grid=list(DISTRIBUTION_TM), cases=rows,
        protocol='Equilibrate max(200tau_c,20tau_m), rounded up to native snapshot cadence; '
                 'sample at least600tau_c. All native sites and equally spaced frames retained. '
                 'Inspect window drift, retain flagged cases, no assumed iid cell/frame samples.',
        comparison={'m': 'Beta on [0,1], matched empirical mean and population variance',
                    'pressure': 'Gaussian on the full real line, matched empirical mean and population variance'},
        output=dict(frame_light=1, ninfo=ninfo, nvideo=0, fields=['pressure', 'm'],
                    steady_frames=sample_steps//ninfo+1))
    immutable_write(root/'manifest.json', json.dumps(manifest, indent=2)+'\n')
    immutable_write(root/'README.md', '# Matched pressure and memory distributions\n\n'
        '96 runs: 4 fixed activities (.3,.5,.7,1) x8 memory times (.3,1,3,6,8,10,20,30) x3 seeds.\n\n'
        'Core physics exactly inherited from mc=.2287, L256 figure. pc=.016838, tau_chi=202.3, '
        'r=.3, tau_c=674.3290333006435 simulation steps. Activity is prescribed (open loop).\n\n'
        'Memory compared with Beta and physical pressure with Gaussian having the same first '
        'and second moments. No pressure recentering, clipped tails, MLE or iid-site inference. '
        'The same simulations supply both fields. All seed outcomes remain visible.\n\n'
        'Warmup max(200tau_c,20tau_m); measurement >=600tau_c; frames every~5tau_c. '
        'Seven native fields retained by frame-light; no video. Analysis: cw_poster_distribution.py.\n')
    return manifest


def pulse_campaign(parent, base, provenance):
    root = parent/PULSE_CAMPAIGN
    rows = []
    tp = round(3*TC)
    wait = round(2000*TC)
    recovery = round(1000*TC)
    for tm in PULSE_TM:
        prep = round(max(100., 10*tm)*TC)
        start = prep+wait
        nsteps = start+tp+recovery
        for init in ('0', '1'):
            for rep in (1, 2, 3):
                stem = f'L256/tm{token(tm)}_init{init}_rep{rep}'
                control = 'controls/'+stem
                for kind in ('control', 'pulse'):
                    name = control if kind == 'control' else 'pulses/tp3/'+stem
                    duration = tp if kind == 'pulse' else 0
                    p = dict(base)
                    p.update({'open-loop': 0, 'seed': 190900+rep, 'chi-seed': 200900+rep,
                              'chi-config': 'uniform', 'chi0': int(init),
                              'm0': .3301 if init == '0' else .0559,
                              'tau-m': tm*TC, 'chi-freeze-steps': prep, 'mem-freeze-steps': 0,
                              'nsteps': nsteps, 'nstart': 0, 'ninfo': round(100*TC),
                              'pmem-pulse-value': .9*PC, 'pmem-pulse-start': start,
                              'pmem-pulse-steps': duration, 'nresponse': 13,
                              'nvideo': 0, 'nvideo-dense': 0, 'video-start': 0,
                              'ntracer': 0, 'tracer-count': 0, 'frame-light': 1})
                    row = write_case(root, name, p, [PULSE_CAMPAIGN,
                        f'Paired {kind}; core mc=.2287, pc=.016838; tau_m/tau_c={tm:g}.',
                        '10% compression-sensing pulse: pc=.016838 -> .0151542 -> .016838.',
                        'The memory threshold changes; no direct mechanical force is added.',
                        f'Frozen-reaction preparation {prep}; full-feedback wait {wait}; pulse {tp} steps.',
                        f'All paired arms end at {nsteps}; >=1000tau_c recovery; exact means every13steps.',
                        'Three signed paired responses are averaged before fitting; retain all fit statuses.'], base)
                    row.update(kind=kind, paired_control=control, tm_over_tc=tm,
                               tp_over_tc=3. if kind == 'pulse' else 0., initialization=init,
                               replicate=rep, seed=190900+rep, chi_seed=200900+rep,
                               L=256, mc=MC, preparation_steps=prep, feedback_wait_steps=wait,
                               pulse_start_steps=start, pulse_duration_steps=duration,
                               recovery_steps_min=recovery, nsteps=nsteps, near_critical=7. <= tm <= 9.)
                    rows.append(row)
    manifest = common_manifest(PULSE_CAMPAIGN, base, provenance)
    manifest['fixed']['pmem_pulse_value'] = .9*PC
    manifest.update(tm_grid=list(PULSE_TM), tp_grid=[3.], initializations=['0', '1'], cases=rows,
        protocol='Both uniform initial histories at every tau_m,3 paired seeds each. '
                 'Prepare max(100tau_c,10tau_m) with chi reaction frozen, then wait2000tau_c '
                 'with feedback. Pulse pc down10% for3tau_c and observe1000tau_c. '
                 'Same-seed controls have identical prehistory and end time. '
                 'Fixed age is not a stationarity guarantee; report drift and residuals, '
                 'never selectively average only successfully fitted seeds.',
        output=dict(nresponse=13, ninfo=round(100*TC), nvideo=0, frame_light=1),
        aggregation='Mean of all3 signed paired pulse-control curves at identical steps before '
                    'relaxation fit. Pointwise seed SEM is descriptive; unresolved remains unresolved.')
    immutable_write(root/'manifest.json', json.dumps(manifest, indent=2)+'\n')
    immutable_write(root/'README.md', '# Dense 10% compression-sensing pulse scan\n\n'
        '300 runs:25 memory times x2 uniform initial histories x3 seeds x(control+pulse).\n\n'
        'g=tau_m/tau_c: '+', '.join(f'{g:g}' for g in PULSE_TM)+'.\n\n'
        'Core physics exactly inherited from the mc=.2287 L256 figure. Baseline pc=.016838; '
        'pulse pc=.0151542 for3tau_c; tau_chi=202.3, r=.3. This changes memory input only.\n\n'
        'Every point uses the same2000tau_c feedback wait and1000tau_c post-pulse observation '
        'after max(100tau_c,10tau_m) frozen preparation. These are finite ages, not imposed '
        'stationarity. Exact global response sampled every13steps and at events; sparse native '
        'fields every~100tau_c; no large video streams.\n\n'
        'Three seeds are averaged as signed pulse-control response before fitting. Each '
        'initialization is separate. All drift/persistent/missing-fit states retained; '
        'no forced zero/infinite relaxation times. Analysis: cw_poster_pulse.py.\n')
    return manifest


def validate(parent, manifests, base):
    total = 0
    for manifest in manifests:
        root = parent/manifest['campaign']
        cases = manifest['cases']
        expected = {root/r['case']/'run.dat' for r in cases}
        if len(expected) != len(cases) or expected != set(root.rglob('*.dat')):
            raise ValueError('Duplicate, missing or stray dat files')
        by_name = {r['case']: r for r in cases}
        groups = Counter()
        for row in cases:
            path = root/row['case']/'run.dat'
            if hashlib.sha256(path.read_bytes()).hexdigest() != row['run_dat_sha256']:
                raise ValueError('Runcard digest mismatch')
            p = parse_card(path.read_text())
            if p != row['expected_parameters'] or any(p[k] != base[k] for k in PHYSICS_KEYS):
                raise ValueError('Core/reference parameter mismatch')
            if not math.isclose(float(p['tau-m'])/TC, row['tm_over_tc'], rel_tol=1e-13):
                raise ValueError('Memory clock mismatch')
            if row.get('kind') == 'pulse':
                q = by_name[row['paired_control']]['expected_parameters']
                if {k for k in p if p[k] != q[k]} != {'pmem-pulse-steps'}:
                    raise ValueError('Paired history differs beyond pulse duration')
                if not math.isclose(1-float(p['pmem-pulse-value'])/PC, .1, abs_tol=1e-14):
                    raise ValueError('Pulse amplitude mismatch')
            key = (row['tm_over_tc'], row.get('activity'), row.get('initialization'), row.get('kind'))
            groups[key] += 1
        if set(groups.values()) != {3}:
            raise ValueError('Every condition must have exactly three seeds')
        total += len(cases)
    if [len(m['cases']) for m in manifests] != [96, 300] or total != 396:
        raise ValueError('Campaign size mismatch')
    return dict(validated_cases=total, source_files=len(manifests[0]['simulator_source_sha256']),
                matched_core_physics=list(PHYSICS_KEYS), distributions=96, pulse=150, control=150)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-root', type=Path, default=REPO/'cases')
    args = parser.parse_args()
    base = parse_card((REPO/BASELINE).read_text())
    provenance = source_provenance()
    manifests = [distribution_campaign(args.out_root, base, provenance),
                 pulse_campaign(args.out_root, base, provenance)]
    print(json.dumps(validate(args.out_root, manifests, base), indent=2))
