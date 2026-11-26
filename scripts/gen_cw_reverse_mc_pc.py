#!/usr/bin/env python3
"""Reversed feedback (m above mc -> PASSIVE): is there bistability anywhere in (mc, pc)?  L128, step law.

20261125/cw_reverse_feedback found one state at every tau_m at the working point (mc = .2287,
pc = .016838). The mean-field reason is the sign of the loop: with f(a) = Prob(P > pc) in the
uniform state of active fraction a, the map a -> Theta(mc - f(a)) is decreasing whenever f grows
with a, and a decreasing map crosses the diagonal once. f(a) was measured in the two uniform states
of that campaign (pos, tau_m = 18, 20 full frames each, L128):

    pc/pc0     -1       0      0.5      1       2            pc0 = .016838
    f(a=1)    .7050   .5280   .4290   .3299   .1575          active:  sigma_P = .0338
    f(a=0)    .9071   .5325   .2431   .0526   .0048          passive: sigma_P = .0121

f grows with a for pc > 0 (the active state has the wider pressure distribution, so it is above a
POSITIVE threshold more often) and FALLS with a for pc < 0 (it is also below a negative threshold
more often). So this scan asks two things:

  pc > 0  (pc/pc0 = 0.5, 1, 2)   mc across and beyond the window [f(0), f(1)]: predicted ONE state
                                 everywhere -- uniform outside the window, mixed with <m> ~ mc inside.
  pc < 0  (pc/pc0 = -1)          the same law, but now the loop is self-reinforcing (the active
                                 state lowers m, and low m keeps the cell active): predicted TWO
                                 branches inside the window at long tau_m, one state at short tau_m.

mc = f_lo + q (f_hi - f_lo) with q = 0.2, 0.5, 0.8 inside the window, q = 1.3 above it, and below it
max(f_lo - 0.3 (f_hi - f_lo), f_lo/2) (the first is negative for pc/pc0 = 1, 2). tau_m = 1 and
14.5 tau_c, both uniform starts (chi0 = 0 active, chi0 = 1 passive; legacy law, chi = 1 PASSIVE),
rep1. m0 is the measured f of the start state at that pc.

Every card is the 20261125/cw_reverse_feedback `rev` card (switch-sign = +1, L128, z0 = .3,
Dbio = .00166667, tau_chi = 202.3 steps, video stride 2, full frame every 100 tau_c) with only
these set: mc, pmem, tau-m, the preparation (chi reaction frozen for max(100, 10 tau_m) tau_c) and
nsteps (+ 2000 tau_c observation), chi0 and m0.

  gen_cw_reverse_mc_pc.py cases/20261126/cw_reverse_mc_pc
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

from gen_cw_reverse_feedback import L, OBS, TC, digest, parameters, replace, token

PC0 = 0.016838
# pc/pc0 -> (f in the uniform active state, f in the uniform passive state), L128, measured
F = {-1: (.7050, .9071), .5: (.4290, .2431), 1: (.3299, .0526), 2: (.1575, .0048)}
Q = ('below', .2, .5, .8, 1.3)
TM = (1, 14.5)
BASE = 'cases/20261125/cw_reverse_feedback/L128/rev_tm1_chi0_rep1/run.dat'
CHANGED = {'mc', 'pmem', 'tau-m', 'chi-freeze-steps', 'nsteps', 'chi0', 'm0'}


def mc_value(r, q):
    lo, hi = sorted(F[r])
    return round(max(lo - .3*(hi - lo), lo/2) if q == 'below' else lo + q*(hi - lo), 4)


def generate(out):
    repo = Path(__file__).resolve().parents[1]
    out = out.resolve()
    campaign = out.relative_to(repo/'cases').as_posix()
    base_body = 'model' + (repo/BASE).read_text().split('\nmodel', 1)[1]
    b = parameters(base_body)
    assert (b['switch-sign'], b['LX'], b['chi-width'], b['pmem-width'], b['seed'], b['chi-seed']) == (
            '1', str(L), '0', '0', '190901', '200901') and 'switch-model' not in b
    rows, pending = [], []
    for r, (fa, fp) in F.items():
        pc = r*PC0
        for q in Q:
            mc = mc_value(r, q)
            assert 0 < mc < 1
            for g in TM:
                prep = round(max(100., 10*g)*TC)
                for start, m0 in (('0', fa), ('1', fp)):
                    text = base_body
                    for key, value in (('mc', f'{mc:g}'), ('pmem', f'{pc:.6g}'), ('tau-m', f'{g*TC:.13g}'),
                                       ('chi-freeze-steps', prep), ('nsteps', prep + OBS),
                                       ('chi0', start), ('m0', f'{m0:g}')):
                        text = replace(text, key, value)
                    state = 'active' if start == '0' else 'passive'
                    where = 'below the window' if q == 'below' else f'q = {q:g}' + (' (above the window)' if q > 1 else '')
                    text = (f'# {campaign}: reversed feedback, L={L}, pc/pc0={r:g}, mc={mc:g} ({where}), '
                            f'tau_m/tau_c={g!r}, start chi={start} ({state}), rep1.\n'
                            '# Legacy step law, chi = 1 PASSIVE; switch-sign = +1, chi* = Theta(m - mc): m above mc -> PASSIVE.\n'
                            f'# pmem = {pc:.6g}; f = Prob(P > pmem) is {fa:g} in the uniform active state and {fp:g} in the passive one '
                            f'(L128, measured), so the window is [{min(fa, fp):g}, {max(fa, fp):g}].\n'
                            f'# Prepare {prep} steps (max(100,10 tau_m) tau_c, chi reaction frozen), observe {OBS} steps (2000 tau_c).\n'
                            + text)
                    a = parameters(text)
                    changed = {key for key in a if a[key] != b.get(key)}
                    assert changed <= CHANGED and set(a) == set(b), changed
                    assert math.isclose(float(a['tau-m'])/TC, g, rel_tol=1e-11)
                    case = f'L{L}/rev_pc{"m" if r < 0 else ""}{token(abs(r))}_mc{token(mc)}_tm{token(g)}_chi{start}_rep1'
                    rows.append({'case': case, 'switch_sign': 1, 'L': L, 'pc_over_pc0': r, 'pmem': float(f'{pc:.6g}'),
                                 'mc': mc, 'q': q, 'f_active': fa, 'f_passive': fp, 'tm_over_tc': g,
                                 'initialization': start, 'start_state': state, 'chi0': float(start), 'm0': m0,
                                 'replicate': 1, 'preparation_steps': prep, 'nsteps': prep + OBS,
                                 'run_dat_sha256': digest(text.encode())})
                    pending.append((out/case/'run.dat', text))
    assert len({r['case'] for r in rows}) == len(rows) == len(F)*len(Q)*len(TM)*2
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-files', 'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {n: digest(subprocess.check_output(['git', 'show', f'HEAD:{n}'], cwd=repo)) for n in names}
    manifest = {
        'campaign': campaign, 'tau_c_reference': TC, 'base_card': BASE,
        'base_card_sha256': digest((repo/BASE).read_bytes()), 'source_base_commit': commit,
        'simulator_source_sha256': hashes, 'pc0': PC0, 'pc_over_pc0': list(F), 'q': list(Q), 'tm_grid': list(TM),
        'f_uniform': {f'{r:g}': {'active': fa, 'passive': fp} for r, (fa, fp) in F.items()},
        'f_source': 'Prob(P > pc) over 20 full frames after release of 20261125/cw_reverse_feedback/L128/'
                    'pos_tm18_chi0_rep1 (uniform active) and pos_tm18_chi1_rep1 (uniform passive, <chi> = .993)',
        'mc': {f'{r:g}': [mc_value(r, q) for q in Q] for r in F},
        'fixed': {'L': L, 'switch_sign': 1, 'zeta0_frac': 0.3, 'Dbio': 0.00166667, 'tau_chi_steps': 202.3,
                  'chi_width': 0, 'pmem_width': 0, 'seed': 190901, 'chi_seed': 200901},
        'convention': 'legacy law: chi=1 PASSIVE, active fraction = 1 - <chi>; chi0=0 active start, chi0=1 passive start; '
                      'm0 = f of the start state at that pc.',
        'protocol': 'Preparation max(100,10 tau_m) tau_c with chi reaction frozen; observe 2000 tau_c; '
                    'terminal window max(500, 20 tau_m) tau_c (reducer).',
        'cases': rows}
    pending.append((out/'manifest.json', json.dumps(manifest, indent=1) + '\n'))
    for path, body in pending:
        if path.exists() and path.read_text() != body:
            raise FileExistsError(f'Refusing to replace differing file: {path}')
    for path, body in pending:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body)
    assert {p.resolve() for p in out.rglob('*.dat')} == {(out/r['case']/'run.dat').resolve() for r in rows}
    print(json.dumps({'campaign': campaign, 'cases': len(rows), 'mc': manifest['mc']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('out', type=Path)
    generate(p.parse_args().out)
