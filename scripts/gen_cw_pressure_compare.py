#!/usr/bin/env python3
"""Generate a fixed-activity L256 pressure-definition comparison campaign."""
import hashlib
import json
from pathlib import Path
import subprocess

REPO = Path(__file__).resolve().parents[1]
CAMPAIGN = '20260915/cw_pressure_compare'
TC = 674.3290333006435


def generate():
    root = REPO / 'cases' / CAMPAIGN
    base = REPO / 'cases/20260904/cw_s3_a/a0p3/run.dat'
    params = dict(line.split('#')[0].strip().split('=', 1)
                  for line in base.read_text().splitlines()
                  if '=' in line.split('#')[0])
    params = {k.strip(): v.strip() for k, v in params.items()}
    params.update(LX='256', LY='256', nsteps='404400', ninfo='6740',
                  nvideo='0', ntracer='0', **{'tracer-count': '0', 'frame-light': '0',
                  'tau-chi': '202.3', 'tau-m': str(10*TC), 'pmem': '0.016838',
                  'mc': '0.2287', 'm0': '0.25'})
    rows = []
    for a in [.3, .45, .6, .65, .8, 1.]:
        for rep in [1, 2, 3]:
            name = f'a{a:g}'.replace('.', 'p') + f'/rep{rep}'
            p = dict(params, seed=str(191500+rep))
            p['zeta-open'] = f'{a*float(p["zeta"]):.12g}'
            body = ('# Fixed-activity pressure comparison; no phenotype feedback into stress.\n'
                    '# p1=pressure; p2=-sigma_bulk; pLB=pressure_lb, all directly serialized.\n'
                    '# Full frames every 6740 steps; analysis window 134800..404400.\n'
                    + ''.join(f'{k:22s} = {v}\n' for k, v in p.items()))
            path = root / name / 'run.dat'
            if path.exists() and path.read_text() != body:
                raise FileExistsError(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(body)
            rows.append(dict(case=name, activity=a, replicate=rep, seed=int(p['seed']),
                             zeta_open=float(p['zeta-open']),
                             run_dat_sha256=hashlib.sha256(body.encode()).hexdigest(),
                             expected_parameters=p))
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=REPO, text=True).splitlines()
    sources = {n: hashlib.sha256(subprocess.check_output(['git', 'show', f'HEAD:{n}'],
                                                        cwd=REPO)).hexdigest() for n in names}
    manifest = dict(campaign=CAMPAIGN, source_base_commit=commit,
                    simulator_source_sha256=sources, cases=rows, tau_c=TC,
                    steady_start=134800, nsteps=404400, ninfo=6740,
                    stationary_gate='Two equal halves: relative changes of u_rms and both '
                    'spatial pressure SDs <=15%; mean changes / pooled SD <=15%. '
                    'All cases remain visible; gate failure is flagged, never hidden.',
                    sigma_definition='Population SD over all sites and stored times in the '
                    'analysis window, computed per seed; trend shows seed mean +/- seed SD.',
                    source_note='Committed simulator only; unrelated local bacterial work excluded.')
    (root / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'{len(rows)} cases in {root}; {404400/TC:.3f} reference tau_c')


if __name__ == '__main__':
    generate()
