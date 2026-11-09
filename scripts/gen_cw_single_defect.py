#!/usr/bin/env python3
"""A single +1/2 in a uniform, non-adaptive layer: the pressure around it under several definitions.

Every card is the 20261022/cw_hill_defect_dense tau_chi = 0.3, tau_m = 8 card with only these set:
  open-loop = 1, zeta-open = zeta * a  (a = 0.3, 0.685, 1; zeta = 0.0666667): uniform activity, chi not fed
                                        back into the stress, so the phenotype dynamics cannot act
  director-config = defect-pair, defect-sep = 128, noise = 0: one +1/2 at x = 64 and one -1/2 at x = 192,
                                        y = 128, no initial director noise (the dynamics has no noise)
  chi-freeze-steps = 0, nstart = 0, nsteps = 10115 (15 tau_c), ninfo = 135 (0.2 tau_c), frame-light = 0:
                                        full frames (with p_LB, sigma_bulk and the stress tensor) from t = 0

  gen_cw_single_defect.py cases/20261109/cw_single_defect
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
SET = (('nsteps', 10115), ('nstart', 0), ('ninfo', 135), ('frame-light', 0), ('chi-freeze-steps', 0),
       ('director-config', 'defect-pair'), ('defect-sep', 128), ('noise', 0))


def digest(b):
    return hashlib.sha256(b).hexdigest()


def replace(text, key, value):
    body, n = re.subn(rf'(?m)^{re.escape(key)}\s*=.*$', f'{key:<21}= {value}', text)
    assert n == 1, key
    return body


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('out', type=Path)
    out = p.parse_args().out
    repo = Path(__file__).resolve().parents[1]; out = out.resolve()
    campaign = out.relative_to(repo/'cases').as_posix()
    body = 'model' + (repo/BASE).read_text().split('\nmodel', 1)[1]
    rows, pending = [], []
    for a in LEVELS:
        text = body
        for k, v in SET:
            text = replace(text, k, v)
        zo = f'{ZETA*a:.7g}'
        text, n = re.subn(r'(?m)^(zeta0-frac\s*=.*\n)', lambda m: m.group(1) + f'{"open-loop":<21}= 1\n{"zeta-open":<21}= {zo}\n', text)
        assert n == 1
        tok = f'{a:g}'.replace('.', 'p')
        text = (f'# {campaign}: single +1/2 (defect pair, sep 128) in a uniform OPEN-LOOP layer, activity a = {a:g} '
                f'(zeta-open = {zo}).\n# Base {BASE}; changed: open-loop, zeta-open, ' + ', '.join(k for k, _ in SET)
                + '.\n# 15 tau_c from t = 0, a full frame every 135 steps (0.2 tau_c).\n' + text)
        case = f'L256/a{tok}'
        rows.append({'case': case, 'activity_a': a, 'zeta_open': float(zo), 'run_dat_sha256': digest(text.encode())})
        pending.append((out/case/'run.dat', text))
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-files', 'src', 'Makefile'], cwd=repo, text=True).splitlines()
    manifest = {'campaign': campaign, 'base_card': BASE, 'base_card_sha256': digest((repo/BASE).read_bytes()),
                'source_base_commit': commit,
                'simulator_source_sha256': {n: digest(subprocess.check_output(['git', 'show', f':{n}'], cwd=repo)) for n in names},
                'changed': dict(SET) | {'open-loop': 1}, 'cases': rows}
    pending.append((out/'manifest.json', json.dumps(manifest, indent=1) + '\n'))
    for path, text in pending:
        if path.exists() and path.read_text() != text:
            raise FileExistsError(path)
    for path, text in pending:
        path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text)
    print(json.dumps({'campaign': campaign, 'cases': [r['case'] for r in rows]}))


if __name__ == '__main__':
    main()
