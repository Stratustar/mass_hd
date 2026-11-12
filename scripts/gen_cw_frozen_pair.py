#!/usr/bin/env python3
"""A +-1/2 pair with Q held fixed, then released (q-freeze-start / q-freeze-steps).

  gen_cw_frozen_pair.py pilot cases/20261112/cw_frozen_pair_pilot

pilot: uniform (open-loop) layer, activity a = 0.3, 0.685, 1 (zeta-open = zeta * a). Every card is the
20261022/cw_hill_defect_dense tau_chi = 0.3, tau_m = 8 card with only these set:
  open-loop = 1, zeta-open:  uniform activity, chi not fed back into the stress
  director-config = defect-pair, defect-sep = 128, noise = 0:  +1/2 at x = 64, -1/2 at x = 192, y = 128
  chi-freeze-steps = 0, nstart = 0, ninfo = 135 (0.2 tau_c), frame-light = 1
  q-freeze-start = 1080 (1.6 tau_c: the cores of the analytic pair relax and the +1/2 starts to move),
  q-freeze-steps = 6750 (10 tau_c with Q held: the flow and the pressure reach the steady state of that Q,
                         the material streams through the pinned +1/2),
  nsteps = 11205: Q released at step 7830 and followed for 5 tau_c.
The pilot measures the steady pressure and velocity around the pinned defects (to place pmem for the
tau_m / tau_chi scan) and checks them against the moving defect after release. m and chi evolve with the
turbulent-layer thresholds and are not used.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

BASE = 'cases/20261022/cw_hill_defect_dense/L256/hn16_wp0p002_tc0p3_tm8_chi1_rep1/run.dat'
ZETA = 0.0666667
LEVELS = (0.3, 0.685, 1.0)
SET = (('nsteps', 11205), ('nstart', 0), ('ninfo', 135), ('frame-light', 1), ('chi-freeze-steps', 0),
       ('director-config', 'defect-pair'), ('defect-sep', 128), ('noise', 0))
ADD = (('q-freeze-start', 1080), ('q-freeze-steps', 6750))


def digest(b):
    return hashlib.sha256(b).hexdigest()


def replace(text, key, value):
    body, n = re.subn(rf'(?m)^{re.escape(key)}\s*=.*$', f'{key:<21}= {value}', text)
    assert n == 1, key
    return body


def pilot(repo, out, campaign):
    body = 'model' + (repo/BASE).read_text().split('\nmodel', 1)[1]
    rows, pending = [], []
    for a in LEVELS:
        text = body
        for k, v in SET:
            text = replace(text, k, v)
        zo = f'{ZETA*a:.7g}'
        text, n = re.subn(r'(?m)^(zeta0-frac\s*=.*\n)', lambda m: m.group(1) + f'{"open-loop":<21}= 1\n{"zeta-open":<21}= {zo}\n', text)
        assert n == 1
        text = text.rstrip('\n') + ''.join(f'\n{k:<21}= {v}' for k, v in ADD) + '\n'
        tok = f'{a:g}'.replace('.', 'p')
        text = (f'# {campaign}: +-1/2 pair (sep 128), uniform OPEN-LOOP layer, activity a = {a:g} (zeta-open = {zo}); '
                f'Q held for steps 1080-7829, then released.\n# Base {BASE}; changed: open-loop, zeta-open, '
                + ', '.join(k for k, _ in SET + ADD) + '.\n# 16.6 tau_c from t = 0, a light frame every 135 steps (0.2 tau_c).\n' + text)
        case = f'L256/a{tok}'
        rows.append({'case': case, 'activity_a': a, 'zeta_open': float(zo), 'run_dat_sha256': digest(text.encode())})
        pending.append((out/case/'run.dat', text))
    return rows, pending, dict(SET + ADD) | {'open-loop': 1}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('stage', choices=['pilot'])
    p.add_argument('out', type=Path)
    a = p.parse_args()
    repo = Path(__file__).resolve().parents[1]; out = a.out.resolve()
    campaign = out.relative_to(repo/'cases').as_posix()
    rows, pending, changed = {'pilot': pilot}[a.stage](repo, out, campaign)
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-files', 'src', 'Makefile'], cwd=repo, text=True).splitlines()
    manifest = {'campaign': campaign, 'stage': a.stage, 'base_card': BASE, 'base_card_sha256': digest((repo/BASE).read_bytes()),
                'source_base_commit': commit,
                'simulator_source_sha256': {n: digest(subprocess.check_output(['git', 'show', f':{n}'], cwd=repo)) for n in names},
                'changed': changed, 'cases': rows}
    pending.append((out/'manifest.json', json.dumps(manifest, indent=1) + '\n'))
    for path, text in pending:
        if path.exists() and path.read_text() != text:
            raise FileExistsError(path)
    for path, text in pending:
        path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text)
    print(json.dumps({'campaign': campaign, 'cases': [r['case'] for r in rows]}))


if __name__ == '__main__':
    main()
