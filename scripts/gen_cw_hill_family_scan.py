#!/usr/bin/env python3
"""Hill-law tau_m scans in two families, 3 seeds, at mc = 0.2287 (L256).

  nscan     Hill exponent n = 4, 8, 16, 32, 64 (wp = 0.002, tau_chi = 0.3 tau_c), plus the step
            limit h-n = 0, wp = 0 as the same-grid reference: does the bifurcation curve approach
            the step-law result as n grows?
  tchiscan  tau_chi = 1, 3, 10, 30 tau_c at n = 16 (tau_chi = 0.3 tau_c is the nscan n = 16
            member and is not duplicated): does it approach the single-cell curve as tau_chi grows?
  nwscan    n and wp together, wp = 0.032/n: (4, .008), (8, .004), (32, .001), (64, .0005).
            With wp held at .002 the nscan converged in n (32 ~ 64) but its passive branch stayed
            at tau_m ~ 10 against 8.5 for the step limit: a fixed wp is a fixed absolute width,
            ~0.17 sigma_P in the passive phase, and it lifts the passive threshold statistic f by
            ~12 %. (16, .002) and the step limit are the nscan members and are not duplicated.
  dscan     diffusion: Dm = Dchi = .00166667/den, den = 3, 10, 30, 100, at n = 16, wp = .002,
            tau_chi = 0.3 tau_c; tau_m = 2, 8, 12, 18, 30 tau_c, 2 seeds, a full frame every
            ~250 tau_c for spatial spectra. den = 1 is run only at tau_m = 30: at 2-18 it is the
            nscan n = 16 member (same seeds; output cadence does not change the dynamics).
            Small D cannot blow up (explicit diffusion has only an upper bound); the question is
            resolution -- the Taylor scale of chi and m is ~2 cells at den = 1 and scales as
            sqrt(D) -- and how the bifurcation moves, since Dm sets the memory width.

Every card is the 20261013/cw_hill_tm_scan card (hill law, beta = 0, Dm = Dchi = .00166667,
z0 = .3, pmem = .016838) with only these set: tau-m, the preparation (chi reaction frozen for
max(100, 10 tau_m) tau_c) and nsteps (+ 2000 tau_c observation), the start, the seeds, h-n, wp,
tau-chi, and the OUTPUT (video-stride 8, a full frame every ~500 tau_c) -- output settings do not
touch the dynamics. Exact full-field means are still written every 337 steps (0.5 tau_c).

  gen_cw_hill_family_scan.py nscan    cases/20261015/cw_hill_nscan
  gen_cw_hill_family_scan.py tchiscan cases/20261015/cw_hill_tchiscan
  gen_cw_hill_family_scan.py nwscan   cases/20261016/cw_hill_nwscan
  gen_cw_hill_family_scan.py dscan    cases/20261016/cw_hill_dscan
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

TC = 674.3290333006435
TM = (0.3, 1, 2, 3.5, 5, 6, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14.5, 16, 18)
OBS = 1348658                          # round(2000 tau_c), as every previous L256 scan
NINFO = 337165                         # full frame every ~500 tau_c (was 67400)
REPS = (1, 2, 3)                       # seed 190900+k, chi-seed 200900+k (rep1 = earlier scans)
STARTS = {'1': '0.3301', '0': '0.0559'}   # chi0 -> m0 (active / passive uniform start)
BASE = 'cases/20261013/cw_hill_tm_scan/L256/mc0p2287_hn16_wp0p002_b0_tm8p684210526_chi1_rep1/run.dat'
FAMILIES = {
    'nscan': [dict(hn=n, wp=.002, tc=.3) for n in (4, 8, 16, 32, 64)] + [dict(hn=0, wp=0, tc=.3)],
    'tchiscan': [dict(hn=16, wp=.002, tc=tc) for tc in (1, 3, 10, 30)],
    'nwscan': [dict(hn=n, wp=.032/n, tc=.3) for n in (4, 8, 32, 64)],
    'dscan': [dict(hn=16, wp=.002, tc=.3, den=1, tm=(30,), reps=(1, 2), ninfo=168582)] +
             [dict(hn=16, wp=.002, tc=.3, den=d, tm=(2, 8, 12, 18, 30), reps=(1, 2), ninfo=168582)
              for d in (3, 10, 30, 100)],
}
DBIO = 0.00166667


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


def generate(family, out):
    repo = Path(__file__).resolve().parents[1]
    out = out.resolve()
    campaign = out.relative_to(repo/'cases').as_posix()
    base = (repo/BASE).read_text()
    base_body = 'model' + base.split('\nmodel', 1)[1]
    b = parameters(base_body)
    assert b['switch-model'] == 'hill' and b['beta'] == '0' and b['Dm'] == b['Dchi'] == '0.00166667'
    rows, pending = [], []
    for mem in FAMILIES[family]:
        member = f'hn{token(mem["hn"])}_wp{token(mem["wp"])}_tc{token(mem["tc"])}'
        den, ninfo = mem.get('den'), mem.get('ninfo', NINFO)
        if den is not None:
            member += f'_D1o{den}'
            dval = f'{DBIO/den:.6g}'
        tau_chi = '202.3' if mem['tc'] == .3 else f'{mem["tc"]*TC:.13g}'
        for g in mem.get('tm', TM):
            prep = round(max(100., 10*g)*TC)
            for start, m0 in STARTS.items():
                for k in mem.get('reps', REPS):
                    text = base_body
                    for key, value in (('tau-m', f'{g*TC:.13g}'), ('chi-freeze-steps', prep),
                                       ('nsteps', prep + OBS), ('chi0', start), ('m0', m0),
                                       ('seed', 190900 + k), ('chi-seed', 200900 + k),
                                       ('h-n', mem['hn']), ('wp', mem['wp']), ('tau-chi', tau_chi),
                                       ('video-stride', 8), ('ninfo', ninfo)) + (
                                       (('Dm', dval), ('Dchi', dval)) if den is not None else ()):
                        text = replace(text, key, value)
                    header = (f'# {campaign}: hill law, {member}, tau_m/tau_c={g!r}, start chi={start}, rep{k}.\n'
                              f'# h-n={mem["hn"]} (0 = step), wp={mem["wp"]} (0 = step), tau_chi={tau_chi} steps; '
                              + (f'beta=0, Dm=Dchi={dval} (= .00166667/{den}), mc=.2287, pmem=.016838, z0=.3.\n'
                                 if den is not None else 'beta=0, Dm=Dchi=.00166667, mc=.2287, pmem=.016838, z0=.3.\n')
                              + f'# Prepare {prep} steps (max(100,10 tau_m) tau_c, chi reaction frozen), observe {OBS} '
                              + f'steps (2000 tau_c). Video stride 8, full frame every ~{round(ninfo/TC)} tau_c.\n')
                    text = header + text
                    a = parameters(text)
                    changed = {key for key in a if a[key] != b.get(key)}
                    assert changed <= {'tau-m', 'chi-freeze-steps', 'nsteps', 'chi0', 'm0', 'seed',
                                       'chi-seed', 'h-n', 'wp', 'tau-chi', 'video-stride', 'ninfo',
                                       'Dm', 'Dchi'}, changed
                    assert math.isclose(float(a['tau-m'])/TC, g, rel_tol=1e-11)
                    case = f'L256/{member}_tm{token(g)}_chi{start}_rep{k}'
                    rows.append({'case': case, 'member': member, 'h_n': mem['hn'], 'wp': mem['wp'],
                                 'tau_chi_over_tc': mem['tc'], 'tau_chi_steps': float(tau_chi),
                                 **({'D_divisor': den, 'Dm': float(dval), 'Dchi': float(dval)} if den is not None else {}),
                                 'L': 256, 'mc': 0.2287, 'tm_over_tc': g, 'initialization': start,
                                 'chi0': float(start), 'replicate': k, 'seed': 190900 + k,
                                 'chi_seed': 200900 + k, 'preparation_steps': prep, 'nsteps': prep + OBS,
                                 'run_dat_sha256': digest(text.encode())})
                    pending.append((out/case/'run.dat', text))
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-files', 'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {n: digest(subprocess.check_output(['git', 'show', f':{n}'], cwd=repo)) for n in names}
    manifest = {
        'campaign': campaign, 'family': family, 'tau_c_reference': TC, 'base_card': BASE,
        'base_card_sha256': digest((repo/BASE).read_bytes()), 'source_base_commit': commit,
        'simulator_source_sha256': hashes, 'tm_grid': list(TM), 'replicates': list(REPS),
        'members': FAMILIES[family], 'initializations': list(STARTS),
        'fixed': {'L': 256, 'mc': 0.2287, 'pmem': 0.016838, 'zeta0_frac': 0.3, 'beta': 0,
                  'Dm': 0.00166667, 'Dchi': 0.00166667},
        'convention': 'switch-model=hill: chi=1 ACTIVE; chi0=1 active start (m0=.3301), chi0=0 passive (m0=.0559).',
        'protocol': 'Preparation max(100,10 tau_m) tau_c with chi reaction frozen; observe 2000 tau_c; '
                    'terminal window max(500, 20 tau_m) tau_c (reducer).',
        'shared_member_note': {'tchiscan': 'tau_chi = 0.3 tau_c at n = 16 is the nscan member hn16_wp0p002_tc0p3',
                               'nwscan': '(n, wp) = (16, .002) and the step limit are the nscan members '
                                         'hn16_wp0p002_tc0p3 and hn0_wp0_tc0p3',
                               'dscan': 'den = 1 at tau_m = 2, 8, 12, 18 is the nscan member hn16_wp0p002_tc0p3 '
                                        '(rep 1-2)'}.get(family, ''),
        'output': {'video_stride': 8, 'nvideo': 337, 'ninfo': NINFO, 'frame_light': 1},
        'cases': rows}
    pending.append((out/'manifest.json', json.dumps(manifest, indent=1) + '\n'))
    for path, body in pending:
        if path.exists() and path.read_text() != body:
            raise FileExistsError(f'Refusing to replace differing file: {path}')
    for path, body in pending:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body)
    assert {p.resolve() for p in out.rglob('*.dat')} == {(out/r['case']/'run.dat').resolve() for r in rows}
    print(json.dumps({'campaign': campaign, 'cases': len(rows),
                      'members': [r for r in sorted({r['member'] for r in rows})]}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('family', choices=tuple(FAMILIES))
    p.add_argument('out', type=Path)
    a = p.parse_args()
    generate(a.family, a.out)
