#!/usr/bin/env python3
"""Finite-size scan of the 50% sensing-threshold pulse around the tau_m transition.

Every card is copied from the L256 cw_poster_pulse50 pulse card at the same tau_m and
initialization; only LX/LY, the seed pair and (for controls) the pulse duration change.
The completed L256 campaign (pulses from cw_poster_pulse50, their reused controls from
cw_poster_pulse10) is listed as the L256 reference and is NOT re-run.

  L128: 12 seeds (noise of a box average ~ 1/L, so 4x the L256 seeds for equal error)
  L512:  3 seeds
  g = tau_m/tau_c: 6 6.5 7 7.25 7.5 7.7 7.9 8.1 8.3 8.5 9 10;  init 0 and 1
  each (L, g, init, seed): one pulse run + one zero-duration control with identical history
"""
import hashlib
import json
from pathlib import Path

from gen_cw_poster_campaigns import (REPO, TC, MC, PC, immutable_write, parse_card,
                                     source_provenance, token)

CAMPAIGN = '20261007/cw_pulse_fss'
SOURCE = '20260917/cw_poster_pulse50'
SCRATCH = '/scratch/helu/mass_hd/cases/'
G = (6., 6.5, 7., 7.25, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 9., 10.)
SIZES = {128: 12, 512: 3}
CHANGED = {'LX', 'LY', 'seed', 'chi-seed', 'pmem-pulse-steps'}


def card_text(values, comments):
    return ''.join(f'# {line}\n' for line in comments) + ''.join(
        f'{key:24s} = {value}\n' for key, value in values.items())


def main():
    root = REPO / 'cases' / CAMPAIGN
    source_manifest = REPO / 'cases' / SOURCE / 'manifest.json'
    cases, reference = [], []
    for g in G:
        for init in '01':
            stem = f'tm{token(g)}_init{init}'
            src = REPO / 'cases' / SOURCE / 'pulses/tp3/L256' / f'{stem}_rep1/run.dat'
            base = parse_card(src.read_text())
            pulse_steps = int(base['pmem-pulse-steps'])
            assert base['LX'] == base['LY'] == '256' and pulse_steps == round(3 * TC)
            assert abs(float(base['tau-m']) - g * TC) < 1e-6 and base['chi0'] == init
            pulse_end = int(base['pmem-pulse-start']) + pulse_steps
            for L, nseed in SIZES.items():
                for rep in range(1, nseed + 1):
                    name = f'{stem}_rep{rep:02d}'
                    for kind in ('pulse', 'control'):
                        p = dict(base)
                        p.update({'LX': str(L), 'LY': str(L), 'seed': str(190900 + rep),
                                  'chi-seed': str(200900 + rep),
                                  'pmem-pulse-steps': str(pulse_steps if kind == 'pulse' else 0)})
                        assert {k for k in p if p[k] != base[k]} <= CHANGED
                        case = f'L{L}/{kind}/{name}'
                        text = card_text(p, [
                            CAMPAIGN + f' -- {kind}, L={L}, tau_m/tau_c={g:g}, init {init}, rep {rep}',
                            f'Copied from {SOURCE}/pulses/tp3/L256/{stem}_rep1; changed only '
                            'LX/LY, seed/chi-seed' + (' and pulse duration (0)' if kind == 'control' else '') + '.',
                            'pc .016838 -> .008419 for 3 tau_c after 100 tau_c frozen prep + 2000 tau_c '
                            'feedback; 1000 tau_c recovery. Pulse/control share the full prehistory.'])
                        immutable_write(root / case / 'run.dat', text)
                        cases.append(dict(case=case, kind=kind, L=L, tm_over_tc=g, initialization=init,
                                          replicate=rep, seed=190900 + rep, chi_seed=200900 + rep,
                                          paired=f'L{L}/{"control" if kind == "pulse" else "pulse"}/{name}',
                                          pulse_end_steps=pulse_end, nsteps=int(p['nsteps']),
                                          run_dat_sha256=hashlib.sha256(text.encode()).hexdigest()))
            for rep in (1, 2, 3):
                for kind, raw in (('pulse', f'{SOURCE}/pulses/tp3/L256/{stem}_rep{rep}'),
                                  ('control', f'{SOURCE}/controls/L256/{stem}_rep{rep}')):
                    reference.append(dict(case=f'L256/{kind}/{stem}_rep{rep:02d}', raw=SCRATCH + raw,
                                          kind=kind, L=256, tm_over_tc=g, initialization=init,
                                          replicate=rep, seed=190900 + rep, chi_seed=200900 + rep,
                                          pulse_end_steps=pulse_end))
    manifest = dict(
        campaign=CAMPAIGN, source_campaign=SOURCE, **source_provenance(),
        source_manifest_sha256=hashlib.sha256(source_manifest.read_bytes()).hexdigest(),
        tau_c=TC, fixed=dict(mc=MC, pmem=PC, pmem_pulse_value=.5 * PC, tau_chi=202.3, r=.3,
                             tp_over_tc=3, feedback_wait_tc=2000, recovery_tc=1000),
        tm_grid=list(G), sizes={str(L): n for L, n in SIZES.items()}, reference_L256_seeds=3,
        cases=cases, reference_L256=reference,
        metrics=dict(
            recovery='T_half: first time the seed-averaged, 5 tau_c-smoothed pulse-minus-control '
                     '<chi> response falls below half its value just after the pulse',
            fluctuation='sigma^2: variance of <chi>(t) in unperturbed controls over the last '
                        '2000 tau_c; S = L^2 sigma^2 for finite-size comparison'))
    immutable_write(root / 'manifest.json', json.dumps(manifest, indent=2) + '\n')
    immutable_write(root / 'README.md', '# Finite-size scan of the 50% pulse response\n\n'
        'L128 (12 seeds) and L512 (3 seeds) at tau_m/tau_c = ' + ', '.join(f'{g:g}' for g in G) +
        ', both uniform initializations; one pulse and one zero-duration control per seed. '
        'L256 is the completed cw_poster_pulse50 campaign (listed in the manifest, not re-run).\n\n'
        'All physics and protocol inherited from cw_poster_pulse50: mc=.2287, pc=.016838, '
        'pulse pc=.008419 for 3 tau_c, r=.3, tau_chi=202.3, 2000 tau_c feedback wait, '
        '1000 tau_c recovery. Only box size and seeds change.\n\n'
        'Metrics (plot/python/confluent_wet/cw_pulse_fss.py): half-recovery time T_half of '
        '<chi> after the pulse, and the stationary fluctuation sigma^2 = Var_t <chi> of the '
        'controls (S = L^2 sigma^2), plus autocorrelation time, Binder cumulant and escape counts.\n')
    print(json.dumps({'campaign': CAMPAIGN, 'new_runs': len(cases),
                      'per_size': {L: sum(c['L'] == L for c in cases) for L in SIZES},
                      'reference_L256_runs': len(reference)}))


if __name__ == '__main__':
    main()
