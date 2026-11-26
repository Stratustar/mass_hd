#!/usr/bin/env python3
"""Interface dynamics: a box prescribed half active (chi = 1), half passive (chi = 0), then released with adaptation on.

One active band of width L/2 centred at x = L/2 and one passive band, two straight boundaries normal to x (periodic box):
    chi(x) = 1/2 [1 - tanh((|x - L/2| - L/4)/w)],  w = 1
m is prescribed with the same step: m_A where chi > 1/2, m_P elsewhere; m_A = .308 and m_P = .064 are the means deep
inside the two bands of the fixed-stripe run 20261124/cw_stripe_toy/lam256_w1 (last 50 tau_c, |d| > 32). Both lie outside
every mc used here, so at t = 0 each band is self-consistent whatever mc is and the three cards differ in mc only.
Q, u and P are the last frame of that fixed-stripe run (100 tau_c of turbulence over this same pattern at z0 = .3), so
the flow already matches the activity step. From t = 0 everything evolves: chi is advected, diffuses and switches.

Three points of the (tau_m, mc) phase diagram (20260909/cw_mc_tm_phase, L256, 2000 tau_c, two uniform starts), all at
tau_m = 20 tau_c so that only mc changes:
    mc = .16   active   both uniform starts end active (the window of history dependence is mc = .20 - .24)
    mc = .22   bistable each uniform start keeps its state
    mc = .28   passive  both uniform starts end passive
Two sets:
    z0p3_*   zeta0-frac = .3, the value the phase diagram was measured at (chi = 0 keeps 30 % of the activity)
    z0_*     zeta0-frac = 0: chi = 0 has no activity at all; the phase diagram does not carry over to this set

--scan tm_tchi (cases/20261126/cw_halfhalf_tm_tchi): the working point mc = .2287, zeta0-frac = .3 kept, and the two times
varied instead, 3 x 3. Picked from the measured <chi>(tau_m, tau_chi) of 20261019/cw_hill_tchifill (two uniform starts):
    tau_m   = 5     one mixed state for every tau_chi (<chi> = .29 - .44, both starts agree)
            = 9     the onset: one state at tau_chi = .3 (<chi> = .71) and 3 (.39), two at tau_chi = 30 (.81 / .13)
            = 14.5  two states for every tau_chi (.96 / .03 - .06)
    tau_chi = .3, 3, 30 tau_c
The start is the same frame (--frame: reuse the init_frame.json already written by the mc scan instead of a new copy).

Base card: cases/20261022/cw_hill_defect_dense/L256/hn16_wp0p002_tc0p3_tm8_chi1_rep1 (hill law, chi = 1 active).
Changed keys: nsteps, nstart, ninfo, chi-freeze-steps, nvideo, tau-m, init-frame, init-frame-phenotype, and mc, zeta0-frac
(--scan mc) or tau-chi (--scan tm_tchi).

  gen_cw_halfhalf.py localtest/20261125/cw_halfhalf [--tau-c 150] [--tau-m 20] [--container-root /work/mass_hd]
  gen_cw_halfhalf.py cases/20261126/cw_halfhalf_tm_tchi --scan tm_tchi --tau-c 300 --container-root /home/helu/mass_hd \
      --frame cases/20261125/cw_halfhalf/L256/init_frame.json
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np

REPO = Path(__file__).resolve().parents[1]
CARD = REPO/'cases/20261022/cw_hill_defect_dense/L256/hn16_wp0p002_tc0p3_tm8_chi1_rep1/run.dat'
START = REPO.parent/'localtest/cases/20261124/cw_stripe_toy/L256/lam256_w1/frame67500.json'
TC = 674.3290333006435
L, W = 256, 1.
M_A, M_P = .308, .064
FLOW = ('QQxx', 'QQyx', 'ux_mat', 'uy_mat', 'pressure')
POINTS = ((.16, 'active'), (.22, 'bistable'), (.28, 'passive'))
Z0S = ((.3, 'z0p3'), (0., 'z0'))
TMS = ((5., 'one mixed state'), (9., 'the onset of the two states'), (14.5, 'two states'))
TCHIS = ((.3, '202.3'), (3., '2022.987099902'), (30., '20229.87099902'))       # tau_c, and steps as the earlier scans wrote them


def replace(text, key, value):
    body, n = re.subn(rf'(?m)^{re.escape(key)}\s*=.*$', f'{key:<21}= {value}', text)
    assert n == 1, key
    return body


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out', type=Path)
    ap.add_argument('--tau-c', type=float, default=150.)
    ap.add_argument('--tau-m', type=float, default=20., help='in units of tau_c')
    ap.add_argument('--ninfo', type=int, default=337, help='steps between light frames (337 = 0.5 tau_c)')
    ap.add_argument('--container-root', default='/work/mass_hd', help='where the repo is mounted when the cards run')
    ap.add_argument('--scan', choices=('mc', 'tm_tchi'), default='mc')
    ap.add_argument('--frame', type=Path, help='an init_frame.json already in the repo to start from (none is written)')
    a = ap.parse_args()
    root = (REPO/a.out if not a.out.is_absolute() else a.out)/'L256'
    root.mkdir(parents=True, exist_ok=True)
    nsteps = int(round(a.tau_c*TC/a.ninfo))*a.ninfo

    src = json.loads(START.read_text())['data']
    x = np.arange(L)
    prof = .5*(1 - np.tanh((np.abs(x - L/2) - L/4)/W))
    chi = np.repeat(prof[:, None], L, axis=1)                # axis 0 of a field is x
    m = np.where(chi > .5, M_A, M_P)
    data = {k: {'type': 'array(double)', 'value': src[k]['value']} for k in FLOW}
    data['chi'] = {'type': 'array(double)', 'value': chi.ravel().tolist()}
    data['m'] = {'type': 'array(double)', 'value': m.ravel().tolist()}
    if a.frame is None:
        (root/'init_frame.json').write_text(json.dumps({'id': 'frame', 'version': 1, 'data': data}))
        frame = root/'init_frame.json'
    else:
        frame = REPO/a.frame
        have = json.loads(frame.read_text())['data']
        assert all(np.array_equal(np.asarray(have[k]['value']), np.asarray(data[k]['value'])) for k in data), 'frame differs'

    text = CARD.read_text()
    base = '\n'.join(ln for ln in text.splitlines() if ln.strip() and not ln.startswith('#'))
    for k, v in (('nsteps', nsteps), ('nstart', 0), ('ninfo', a.ninfo), ('chi-freeze-steps', 0), ('nvideo', 0)):
        base = replace(base, k, v)
    extra = {'init-frame': f'{a.container_root}/{frame.relative_to(REPO)}', 'init-frame-phenotype': 1}
    common = [f'# Start: Q, u, P of the last frame of localtest/20261124/cw_stripe_toy/L256/lam256_w1 (turbulence over this same pattern held fixed); chi and m prescribed, boundaries at x = {L//4} and {3*L//4}, tanh scale {W:g}.',
              f'# {nsteps} steps = {nsteps/TC:.1f} tau_c, a light frame every {a.ninfo} steps ({a.ninfo/TC:.2f} tau_c).']
    if a.scan == 'tm_tchi':
        manifest = dict(study=str(a.out), scan='tm_tchi', nsteps=nsteps, ninfo=a.ninfo, tau_c=TC, m_active=M_A, m_passive=M_P, tanh_scale=W,
                        start=str(START.relative_to(REPO.parent)), frame=str(frame.relative_to(REPO)), chi_profile=prof.tolist(), cases={})
        for tchi, tchi_steps in TCHIS:
            for tm, what in TMS:
                name = f'tc{tchi:g}_tm{tm:g}'.replace('.', 'p')
                d = root/name
                d.mkdir(parents=True, exist_ok=True)
                body = replace(replace(base, 'tau-m', f'{tm*TC:.9f}'), 'tau-chi', tchi_steps)
                head = [f'# 20261126/cw_halfhalf_tm_tchi/{name}: half the box prescribed active (chi = 1, m = {M_A}), half passive (chi = 0, m = {M_P}), then released with adaptation on.',
                        f'# tau_m = {tm:g} tau_c, tau_chi = {tchi:g} tau_c at the working point (mc = .2287, zeta0-frac = .3): {what} in the uniform-start scan 20261019/cw_hill_tchifill.',
                        common[0],
                        f'# Base: {CARD.relative_to(REPO)}; changed: nsteps, nstart, ninfo, chi-freeze-steps, nvideo, tau-m, tau-chi, init-frame, init-frame-phenotype.',
                        common[1]]
                (d/'run.dat').write_text('\n'.join(head) + '\n' + body + '\n' + '\n'.join(f'{k:<21}= {v}' for k, v in extra.items()) + '\n')
                manifest['cases'][name] = dict(tau_m=tm, tau_chi=tchi, regime=what)
                print(f'{name:14s} tau_m {tm:g}, tau_chi {tchi:g}: {what}')
        (root/'manifest.json').write_text(json.dumps(manifest, indent=1))
        print(f'{len(manifest["cases"])} cases under {root}; {nsteps} steps = {nsteps/TC:.1f} tau_c each')
        return
    base = replace(base, 'tau-m', f'{a.tau_m*TC:.9f}')
    manifest = dict(study=str(a.out), nsteps=nsteps, ninfo=a.ninfo, tau_c=TC, tau_m_over_tc=a.tau_m, m_active=M_A, m_passive=M_P,
                    tanh_scale=W, start=str(START.relative_to(REPO.parent)), chi_profile=prof.tolist(), cases={})
    tm = f'{a.tau_m:g}'.replace('.', 'p')
    for z0, ztag in Z0S:
        for mc, phase in POINTS:
            name = f'{ztag}_tm{tm}_mc{mc:g}'.replace('0.', '0p')
            d = root/name
            d.mkdir(parents=True, exist_ok=True)
            body = replace(replace(base, 'mc', mc), 'zeta0-frac', f'{z0:g}')
            where = (f'{phase} region of the (tau_m, mc) phase diagram' if z0 == .3 else
                     f'mc of the {phase} point; the phase diagram was measured at zeta0-frac = .3 and does not carry over')
            head = [f'# 20261125/cw_halfhalf/{name}: half the box prescribed active (chi = 1, m = {M_A}), half passive (chi = 0, m = {M_P}), then released with adaptation on.',
                    f'# tau_m = {a.tau_m:g} tau_c, mc = {mc:g}, zeta0-frac = {z0:g}: {where}.',
                    f'# Start: Q, u, P of the last frame of localtest/20261124/cw_stripe_toy/L256/lam256_w1 (turbulence over this same pattern held fixed); chi and m prescribed, boundaries at x = {L//4} and {3*L//4}, tanh scale {W:g}.',
                    f'# Base: {CARD.relative_to(REPO)}; changed: nsteps, nstart, ninfo, chi-freeze-steps, nvideo, tau-m, mc, zeta0-frac, init-frame, init-frame-phenotype.',
                    f'# {nsteps} steps = {nsteps/TC:.1f} tau_c, a light frame every {a.ninfo} steps ({a.ninfo/TC:.2f} tau_c).']
            (d/'run.dat').write_text('\n'.join(head) + '\n' + body + '\n' + '\n'.join(f'{k:<21}= {v}' for k, v in extra.items()) + '\n')
            manifest['cases'][name] = dict(mc=mc, z0=z0, phase=phase)
            print(f'{name:22s} mc {mc:g}, z0 {z0:g}, {phase}')
    (root/'manifest.json').write_text(json.dumps(manifest, indent=1))
    print(f'{len(manifest["cases"])} cases under {root}; {nsteps} steps = {nsteps/TC:.1f} tau_c each')


if __name__ == '__main__':
    main()
