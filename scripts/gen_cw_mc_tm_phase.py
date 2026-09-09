#!/usr/bin/env python3
"""Generate the long-time, two-initialization empirical phase-map campaign."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

TC = 674.3290333006435
TM = [0.3, 1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 30, 50]
MC = [0.06, 0.10, 0.14, 0.16, 0.18, 0.20, 0.21, 0.22, 0.24,
      0.26, 0.28, 0.30, 0.34, 0.38]
CONTROL_TM = [6, 12, 20, 30, 50]
CONTROL_MC = [0.18, 0.21, 0.24]


def token(x):
    return f'{x:g}'.replace('.', 'p')


def generate(root):
    repo = Path(__file__).resolve().parents[1]
    base = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    source_hashes = {name: hashlib.sha256(subprocess.check_output(
        ['git', 'show', f'HEAD:{name}'], cwd=repo)).hexdigest() for name in names}
    rows = []
    for L, tms, mcs in [(256, TM, MC), (384, CONTROL_TM, CONTROL_MC)]:
        for tm in tms:
            prep = round(max(100., 10*tm)*TC)
            observation = round(max(2000., 60*tm)*TC)
            for mc in mcs:
                for chi in (0, 1):
                    for rep in (1, 2):
                        seed = 190900 + rep + (1000 if L == 384 else 0)
                        rel = f'L{L}/mc{token(mc)}_tm{token(tm)}_chi{chi}_rep{rep}'
                        params = {
                            'model': 'confluent-wet', 'nsteps': prep+observation,
                            'nstart': 0, 'ninfo': 67400, 'LX': L, 'LY': L, 'bc': 0,
                            'seed': seed, 'Gamma': 0.05, 'xi': 0.4, 'CC': 0.00208333,
                            'LL': 0.01666667, 'tau': 2.0, 'rho': 40, 'friction': 0.0,
                            'zeta': 0.06666667, 'zeta0-frac': 0.3,
                            'angle': 0, 'noise': 0.05, 'initial-order': 1,
                            'director-config': 'uniform', 'defect-sep': 0,
                            'Dbio': 0.00166667, 'chi-config': 'uniform', 'chi0': chi,
                            'chi-seed': seed+10000, 'tau-chi': 202.3,
                            'chi-width': 0.0, 'switch-sign': -1, 'mc': mc,
                            'm0': 0.3301 if chi == 0 else 0.0559,
                            'chi-freeze-steps': prep, 'mem-freeze-steps': 0,
                            'tau-m': tm*TC, 'pmem': 0.016838, 'pmem-width': 0.0,
                            'ntracer': 0, 'tracer-count': 0,
                            'nvideo': 337, 'video-stride': 8,
                            'video-p-scale': 0.20206, 'video-u-scale': 0.04538,
                            'frame-light': 1,
                        }
                        header = (
                            '# 20260909 cw_mc_tm_phase: finite-size, finite-time phase map.\n'
                            f'# tau_c reference = {TC:.13f} steps, inherited from L=500.\n'
                            '# Same physical coefficients and pmem=0.016838 as cw_s3_chi_scan.\n'
                            '# Prepare uniform chi with reaction frozen; flow and m evolve.\n'
                            f'# Release chi after {prep} steps; then observe {observation} steps.\n'
                            '# Preparation >=100 tau_c and >=10 tau_m; NOT counted in observation.\n'
                            '# Observation >=2000 tau_c and >=60 tau_m.\n'
                            '# Tail window >=500 tau_c and >=20 tau_m.\n'
                            '# Two independent Q seeds; a seed is shared across parameters and\n'
                            '# chi arms within its replicate (paired initial perturbations).\n'
                            '# Exact full-field averages are in video_meta.csv; pixels are coarse.\n'
                        )
                        body = header + ''.join(f'{k:<20} = {v:.13g}\n' if isinstance(v, float)
                                               else f'{k:<20} = {v}\n' for k, v in params.items())
                        path = root/rel/'run.dat'
                        if path.exists() and path.read_text() != body:
                            raise FileExistsError(f'Refusing to overwrite different runcard: {path}')
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text(body)
                        rows.append({'case': rel, 'L': L, 'tm_over_tc': tm, 'mc': mc,
                                     'chi0': chi, 'replicate': rep, 'seed': seed,
                                     'preparation_steps': prep, 'observation_steps': observation,
                                     'nsteps': prep+observation,
                                     'run_dat_sha256': hashlib.sha256(body.encode()).hexdigest()})
    manifest = {'campaign': '20260909/cw_mc_tm_phase', 'tau_c_reference': TC,
                'source_base_commit': base, 'simulator_source_sha256': source_hashes,
                'fixed': {'pmem': 0.016838, 'tau_chi': 202.3, 'r': 0.3},
                'grids': {'L256': {'tm': TM, 'mc': MC},
                          'L384': {'tm': CONTROL_TM, 'mc': CONTROL_MC}},
                'replicates': 2, 'initializations': [0, 1], 'cases': rows,
                'protocol': 'Uniform chi held for max(100 tau_c,10 tau_m) while flow and m '
                'equilibrate, then feedback released for max(2000 tau_c,60 tau_m). '
                'Reference tau_c and physical pmem are frozen across sizes; no L-specific retuning.',
                'classification': {'tail_tc_min': 500, 'tail_tau_m_min': 20,
                    'drift_tolerance': 0.03, 'four_block_range_tolerance': 0.06,
                    'initialization_agreement_tolerance': 0.05,
                    'replicate_range_tolerance': 0.08,
                    'paired_bistable_separation_min': 0.20,
                    'between_arm_gap_min': 0.10},
                'interpretation': 'Finite-time persistence of initialization dependence, '
                'not proof of infinite-time or thermodynamic bistability.'}
    assert len(rows) == 788
    assert len({r['case'] for r in rows}) == len(rows)
    (root/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (root/'README.md').write_text(
        '# Long-time empirical mc / tau_m phase map\n\n'
        'L=256: 13 memory times × 14 thresholds × 2 chi starts × 2 seeds = 728 runs.\n'
        'L=384: 5 memory times × 3 thresholds × 2 chi starts × 2 seeds = 60 controls.\n\n'
        'See manifest.json for exact grids, fixed coefficients, source hashes, timing, '
        'seeds and diagnostic tolerances. Each run prepares its own matching flow and '
        'memory with uniform chi frozen, then releases feedback. Preparation time is '
        'excluded from observation and terminal statistics. There is no simulation code change.\n\n'
        'Per-run analysis: plot/python/confluent_wet/cw_mc_tm_scan.py via PLOT_SCRIPT. '
        'Aggregate using --summary and --manifest. The summary preserves missing, incomplete, '
        'drifting, switching and seed-sensitive cases instead of forcing a phase label.\n')
    print(json.dumps({'cases': len(rows), 'L256': 728, 'L384': 60,
                      'source_files': len(source_hashes), 'manifest': str(root/'manifest.json')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out', type=Path)
    generate(parser.parse_args().out)
