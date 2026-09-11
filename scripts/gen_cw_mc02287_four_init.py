#!/usr/bin/env python3
"""Generate the mc=0.2287, four-initialization continuation of the phase scan."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess

TC = 674.3290333006435
MC = 0.2287
TM = [round(float(Fraction(3, 10) + i * Fraction(177, 190)), 9) for i in range(20)]
STARTS = ('0', '1', 'stripe01', 'noise')


def token(x):
    return f'{x:.9f}'.rstrip('0').rstrip('.').replace('.', 'p')


def generate(out, noise_std=.1, noise_length=7.02):
    if not (0 < noise_std <= .5 and 0 <= noise_length <= 256):
        raise ValueError('Invalid noise standard deviation or smoothing length')
    repo = Path(__file__).resolve().parents[1]
    baseline = repo/'cases/20260909/cw_mc_tm_phase/L256/mc0p21_tm20_chi0_rep1/run.dat'
    base = {}
    for line in baseline.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            key, value = line.split('=', 1)
            base[key.strip()] = value.strip()
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {n: hashlib.sha256(subprocess.check_output(['git', 'show', f'HEAD:{n}'], cwd=repo)).hexdigest()
              for n in names}
    rows = []
    for tm in TM:
        prep = round(max(100., 10*tm)*TC)
        observation = round(max(2000., 60*tm)*TC)
        for start in STARTS:
            for rep in (1, 2):
                seed = 190900+rep  # Same paired Q seeds as the L256 phase scan.
                params = dict(base)
                chi0 = float(start) if start in ('0', '1') else .5
                m0 = .3301 if start == '0' else .0559 if start == '1' else .2503
                params.update({'nsteps': prep+observation, 'mc': MC, 'tau-m': tm*TC,
                    'chi-freeze-steps': prep, 'seed': seed, 'chi-seed': seed+10000,
                    'chi-config': 'stripe' if start == 'stripe01' else 'noise' if start == 'noise' else 'uniform',
                    'chi0': chi0, 'm0': m0, 'chi-noise': noise_std if start == 'noise' else 0.,
                    'chi-length': noise_length if start == 'noise' else 0.})
                if start == 'stripe01':
                    # In the existing implementation the LEFT side uses *_hi.
                    # Deliberately assign that slot chi=0 to make left=0, right=1.
                    params.update({'chi-hi': 0, 'chi-lo': 1, 'm-hi': .3301, 'm-lo': .0559})
                name = f'L256/mc0p2287_tm{token(tm)}_chi{start}_rep{rep}'
                header = (
                    '# 20260911 cw_mc02287_four_init; same L256 baseline and timing as phase scan.\n'
                    f'# tau_m/tau_c={tm:.9f}; tau_c={TC:.13f} steps; mc={MC}.\n'
                    f'# Initial condition at t=0: {start}; paired realization {rep}.\n'
                    f'# Prepare {prep} steps with chi REACTION off; then observe {observation} steps.\n'
                    '# Advection/diffusion remain on during preparation; mixed patterns evolve.\n'
                    '# The release distribution is measured; do not equate it to the t=0 pattern.\n'
                    '# Stripe: left chi=0,m=.3301; right chi=1,m=.0559 (two periodic interfaces).\n'
                    f'# Noise: correlated Gaussian, mean .5, std {noise_std}, length {noise_length}; clamp to [0,1].\n'
                    '# Noise m starts uniform at .2503, the old a=.65 duty estimate; preparation equilibrates it.\n'
                    '# Exact full-field statistics every 337 steps; final window 500 tau_c.\n')
                body = header + ''.join(f'{k:<20} = {v:.13g}\n' if isinstance(v, float)
                                        else f'{k:<20} = {v}\n' for k, v in params.items())
                path = out/name/'run.dat'
                if path.exists() and path.read_text() != body:
                    raise FileExistsError(f'Refusing to overwrite a different runcard: {path}')
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(body)
                rows.append({'case': name, 'L': 256, 'mc': MC, 'tm_over_tc': tm,
                    'initialization': start, 'chi0': chi0, 'replicate': rep, 'seed': seed,
                    'chi_seed': seed+10000, 'preparation_steps': prep,
                    'observation_steps': observation, 'nsteps': prep+observation,
                    'run_dat_sha256': hashlib.sha256(body.encode()).hexdigest()})
    manifest = {'campaign': '20260911/cw_mc02287_four_init', 'tau_c_reference': TC,
        'baseline_runcard': str(baseline.relative_to(repo)), 'source_base_commit': commit,
        'simulator_source_sha256': hashes,
        'fixed': {'L': 256, 'mc': MC, 'pmem': .016838, 'tau_chi': 202.3, 'r': .3},
        'tm_grid': TM, 'tm_grid_rule': '20 linearly spaced points from 0.3 to 18, rounded to 9 decimals',
        'initializations': list(STARTS), 'replicates': 2,
        'noise': {'distribution': 'correlated Gaussian, centered and normalized before clipping to[0,1]',
                  'mean': .5, 'std': noise_std, 'smoothing_length_lattice': noise_length, 'm0': .2503},
        'stripe': {'left': {'chi': 0, 'm': .3301}, 'right': {'chi': 1, 'm': .0559}},
        'protocol': 'Same phase-scan preparation max(100tau_c,10tau_m) with chi reaction frozen, '
                    'flow/memory/advection/diffusion running. Mixed configurations describe t=0, '
                    'not release. Observe max(2000tau_c,60tau_m); terminal window max(500tau_c,20tau_m).',
        'statistics': {'observable': 'exact full-domain mean chi', 'tail_tc_min': 500,
                       'tail_tau_m_min': 20, 'all_start_agreement_tolerance': .05,
                       'replicate_range_tolerance': .08, 'paired_start_separation_min': .20,
                       'between_start_range_gap_min': .10},
        'interpretation': 'Finite-time outcomes of four initial preparations; no assumption that equal effective-potential depths imply spatial coexistence.',
        'cases': rows}
    assert len(TM) == 20 and TM[0] == .3 and TM[-1] == 18
    assert len(rows) == 160 and len({r['case'] for r in rows}) == 160
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (out/'README.md').write_text(
        '# mc=0.2287: four-initialization memory scan\n\n'
        'L=256, 20 linear memory times from 0.3 to 18 reference tau_c, four starts, two paired seeds: 160 runs.\n\n'
        'Starts at t=0: all 0, all 1, left 0/right 1 stripe, correlated Gaussian noise around 0.5. '
        f'Noise std={noise_std}, smoothing length={noise_length} lattice units. '
        'Noise m starts uniform at 0.2503; stripe memory matches its two initial activity branches.\n\n'
        'The phase-scan preparation protocol is retained. Only the chi reaction is frozen; '
        'advection/diffusion can mix stripe/noise before feedback release. The analysis records '
        'the actual release mean and spatial standard deviation.\n\n'
        'After preparation, every run observes 2000 tau_c (1,348,658 steps). '
        'Preparation lasts 100–180 tau_c; total 1,416,091–1,470,037 steps. Tail: 500 tau_c. '
        'All physical coefficients and output cadence are inherited from the previous L256 phase scan.\n\n'
        'Per-run analysis: cw_four_init_scan.py; aggregate with --summary --manifest. '
        'Each initialization remains separately labeled; noise/stripe are never folded into chi0.\n')
    print(json.dumps({'cases': len(rows), 'tm': TM, 'nsteps_range': [min(r['nsteps'] for r in rows),
                      max(r['nsteps'] for r in rows)], 'manifest': str(out/'manifest.json')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out', type=Path)
    parser.add_argument('--noise-std', type=float, default=.1)
    parser.add_argument('--noise-length', type=float, default=7.02)
    args = parser.parse_args()
    generate(args.out, args.noise_std, args.noise_length)
