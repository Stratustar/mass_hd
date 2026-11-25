#!/usr/bin/env python3
"""Reversed feedback (compression suppresses the activity) against the usual loop, tau_m scan, L128.

The usual loop is self-reinforcing: pressure above pmem writes the memory m, m above mc makes the
cell ACTIVE, and the active turbulence is what produces pressure above pmem. Here the middle arrow
is reversed -- m above mc makes the cell PASSIVE -- and nothing else is touched. The question is
whether the tau_m transition (one state below ~8 tau_c, two branches above) survives.

No new executable is needed. The hill law fixes the sign, but the legacy step law carries it as
switch-sign, with chi = 1 PASSIVE (zeta_eff = zeta [z0 + (1 - z0)(1 - chi)]):

  pos   switch-sign = -1   chi* = Theta(mc - m)   high m -> chi = 0 -> active    (the usual loop)
  rev   switch-sign = +1   chi* = Theta(m - mc)   high m -> chi = 1 -> passive   (reversed)

`pos` is the law of the 2026-09 campaigns and the exact mirror of the hill step limit
(h-n = 0, wp = 0); it is re-run here as the same-size control, because at L128 the branches are
noisier and escapes more frequent than at L256. The two arms differ in switch-sign ALONE.

Every card is the 20260911/cw_mc02287_four_init card (legacy step law, mc = .2287, pmem = .016838,
z0 = .3, Dbio = .00166667, tau_chi = 202.3 steps = 0.3 tau_c) with only these set: LX/LY = 128,
switch-sign, tau-m, the preparation (chi reaction frozen for max(100, 10 tau_m) tau_c) and nsteps
(+ 2000 tau_c observation), the start, the seeds, and the output (video-stride 2 -> 64 x 64, a full
frame every 100 tau_c). Exact full-field means are written every 337 steps (0.5 tau_c).

Starts are the two uniform PHYSICAL states in both arms: chi0 = 0 active (m0 = .3301) and
chi0 = 1 passive (m0 = .0559). In `rev` neither is self-consistent (the active state writes
m > mc, which asks for passive, and vice versa), so both leave at release.

  gen_cw_reverse_feedback.py cases/20261125/cw_reverse_feedback
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

TC = 674.3290333006435
L = 128
TM = (0.3, 1, 2, 3.5, 5, 6, 7, 8, 9, 10, 12, 14.5, 18)
OBS = 1348658                              # round(2000 tau_c), as every L256 scan
NINFO = 67433                              # full frame every 100 tau_c
REPS = (1, 2)                              # seed 190900+k, chi-seed 200900+k
STARTS = {'0': '0.3301', '1': '0.0559'}    # legacy chi0 -> m0: active / passive uniform start
ARMS = {'pos': -1, 'rev': 1}               # arm -> switch-sign
BASE = 'cases/20260911/cw_mc02287_four_init/L256/mc0p2287_tm0p3_chi0_rep1/run.dat'
CHANGED = {'LX', 'LY', 'switch-sign', 'tau-m', 'chi-freeze-steps', 'nsteps', 'chi0', 'm0', 'seed',
           'chi-seed', 'video-stride', 'ninfo'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def parameters(text):
    return dict(tuple(x.strip() for x in line.split('=', 1))
                for line in text.splitlines() if line.strip() and not line.lstrip().startswith('#'))


def token(x):
    return f'{x:g}'.replace('.', 'p')


def replace(text, key, value):
    body, n = re.subn(rf'(?m)^{re.escape(key)}\s*=.*$', f'{key:<21}= {value}', text)
    assert n == 1, key
    return body


def generate(out):
    repo = Path(__file__).resolve().parents[1]
    out = out.resolve()
    campaign = out.relative_to(repo/'cases').as_posix()
    base_body = 'model' + (repo/BASE).read_text().split('\nmodel', 1)[1]
    b = parameters(base_body)
    assert (b['switch-sign'], b['chi-width'], b['pmem-width'], b['mc'], b['pmem'], b['zeta0-frac'],
            b['Dbio'], b['tau-chi'], b['nvideo']) == (
            '-1', '0', '0', '0.2287', '0.016838', '0.3', '0.00166667', '202.3', '337')
    assert 'switch-model' not in b          # legacy law
    rows, pending = [], []
    for arm, sign in ARMS.items():
        for g in TM:
            prep = round(max(100., 10*g)*TC)
            for start, m0 in STARTS.items():
                for k in REPS:
                    text = base_body
                    for key, value in (('LX', L), ('LY', L), ('switch-sign', sign),
                                       ('tau-m', f'{g*TC:.13g}'), ('chi-freeze-steps', prep),
                                       ('nsteps', prep + OBS), ('chi0', start), ('m0', m0),
                                       ('seed', 190900 + k), ('chi-seed', 200900 + k),
                                       ('video-stride', 2), ('ninfo', NINFO)):
                        text = replace(text, key, value)
                    law = ('chi* = Theta(m - mc): compression -> PASSIVE (reversed feedback)' if sign == 1
                           else 'chi* = Theta(mc - m): compression -> ACTIVE (the usual loop, control)')
                    state = 'active' if start == '0' else 'passive'
                    text = (f'# {campaign}: {arm}, L={L}, tau_m/tau_c={g!r}, start chi={start} ({state}), rep{k}.\n'
                            f'# Legacy step law, chi = 1 PASSIVE; switch-sign = {sign:+d}, {law}.\n'
                            '# mc=.2287, pmem=.016838, z0=.3, Dbio=.00166667, tau_chi=202.3 steps. The arms differ in switch-sign alone.\n'
                            f'# Prepare {prep} steps (max(100,10 tau_m) tau_c, chi reaction frozen), observe {OBS} steps (2000 tau_c).\n'
                            '# Full-field means every 337 steps; video stride 2 (64 x 64); full frame every 100 tau_c.\n'
                            + text)
                    a = parameters(text)
                    changed = {key for key in a if a[key] != b.get(key)}
                    assert changed <= CHANGED and set(a) == set(b), changed
                    assert math.isclose(float(a['tau-m'])/TC, g, rel_tol=1e-11)
                    case = f'L{L}/{arm}_tm{token(g)}_chi{start}_rep{k}'
                    rows.append({'case': case, 'arm': arm, 'switch_sign': sign, 'L': L, 'mc': 0.2287,
                                 'tm_over_tc': g, 'initialization': start, 'start_state': state,
                                 'chi0': float(start), 'm0': float(m0), 'replicate': k,
                                 'seed': 190900 + k, 'chi_seed': 200900 + k, 'preparation_steps': prep,
                                 'nsteps': prep + OBS, 'run_dat_sha256': digest(text.encode())})
                    pending.append((out/case/'run.dat', text))
    # The two arms are the same cards up to switch-sign.
    for (pa, ta), (pb, tb) in zip(pending[:len(pending)//2], pending[len(pending)//2:]):
        da, db = parameters(ta), parameters(tb)
        assert {key for key in da if da[key] != db[key]} == {'switch-sign'}, (pa, pb)
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-files', 'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {n: digest(subprocess.check_output(['git', 'show', f'HEAD:{n}'], cwd=repo)) for n in names}
    manifest = {
        'campaign': campaign, 'tau_c_reference': TC, 'base_card': BASE,
        'base_card_sha256': digest((repo/BASE).read_bytes()), 'source_base_commit': commit,
        'simulator_source_sha256': hashes, 'tm_grid': list(TM), 'replicates': list(REPS),
        'arms': {'pos': {'switch_sign': -1, 'law': 'chi* = Theta(mc - m), compression -> active (usual loop)'},
                 'rev': {'switch_sign': 1, 'law': 'chi* = Theta(m - mc), compression -> passive (reversed)'}},
        'initializations': list(STARTS),
        'fixed': {'L': L, 'mc': 0.2287, 'pmem': 0.016838, 'zeta0_frac': 0.3, 'Dbio': 0.00166667,
                  'tau_chi_steps': 202.3, 'chi_width': 0, 'pmem_width': 0},
        'convention': 'legacy law: chi=1 PASSIVE, active fraction = 1 - <chi>; chi0=0 active start '
                      '(m0=.3301), chi0=1 passive start (m0=.0559), in both arms.',
        'protocol': 'Preparation max(100,10 tau_m) tau_c with chi reaction frozen; observe 2000 tau_c; '
                    'terminal window max(500, 20 tau_m) tau_c (reducer).',
        'output': {'video_stride': 2, 'nvideo': 337, 'ninfo': NINFO, 'frame_light': 1},
        'cases': rows}
    pending.append((out/'manifest.json', json.dumps(manifest, indent=1) + '\n'))
    for path, body in pending:
        if path.exists() and path.read_text() != body:
            raise FileExistsError(f'Refusing to replace differing file: {path}')
    for path, body in pending:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body)
    assert {p.resolve() for p in out.rglob('*.dat')} == {(out/r['case']/'run.dat').resolve() for r in rows}
    print(json.dumps({'campaign': campaign, 'cases': len(rows), 'arms': list(ARMS), 'tm': list(TM)}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('out', type=Path)
    generate(p.parse_args().out)
