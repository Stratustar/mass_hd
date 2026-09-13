#!/usr/bin/env python3
"""Complete the 160 held pulses at a prescribed age, without stationarity gates."""
import hashlib
import json
from pathlib import Path
import subprocess

TC = 674.3290333006435
CAMPAIGN = '20260913/cw_pmem_pulse_long_wait'


def main():
    repo = Path(__file__).resolve().parents[1]
    old_root = repo / 'cases/20260913/cw_pmem_pulse'
    old = json.loads((old_root / 'manifest.json').read_text())
    root = repo / 'cases' / CAMPAIGN
    rows = []
    for item in old['cases']:
        near = item['near_critical']
        if not near and not (item['tm_over_tc'] in (12., 20.) and item['initialization'] == '1'):
            continue
        row = dict(item)
        text = (old_root / item['case'] / 'run.dat').read_text()
        if near:
            p = {k.strip(): v.strip() for line in text.splitlines()
                 if line.strip() and not line.startswith('#') for k, v in [line.split('=', 1)]}
            prep = item['preparation_steps']
            wait = round(2000 * TC)
            start = prep + wait
            nsteps = start + round(10 * TC) + round(500 * TC)
            p.update({'pmem-pulse-start': start, 'nsteps': nsteps,
                      'video-start': start - round(20 * TC),
                      'video-dense-start': start - round(20 * TC),
                      'video-dense-end': start + round(110 * TC)})
            text = (f'# {CAMPAIGN}: fixed-age pulse, no drift/velocity submission veto.\n'
                    '# Frozen preparation 100 tau_c; feedback-on wait 2000 tau_c; recovery >=500 tau_c.\n'
                    '# Same-seed controls share the entire prehistory. pmem .016838 -> 0 -> .016838.\n'
                    + ''.join(f'{k:<24} = {v}\n' for k, v in p.items()))
            row.update(feedback_wait_steps=wait, pulse_start_steps=start, nsteps=nsteps)
        elif item['kind'] == 'control':
            row['reuse_from'] = old['campaign'] + '/' + item['case']
        row['run_dat_sha256'] = hashlib.sha256(text.encode()).hexdigest()
        path = root / row['case'] / 'run.dat'
        if path.exists() and path.read_text() != text:
            raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        rows.append(row)
    assert len(rows) == 200
    assert sum(r['kind'] == 'pulse' for r in rows) == 160
    assert sum('reuse_from' in r for r in rows) == 8
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    manifest = {k: old[k] for k in ('fixed', 'simulator_source_sha256', 'tp_grid', 'replicates', 'output')}
    manifest.update(campaign=CAMPAIGN, source_base_commit=commit,
        tm_grid=[7.7, 7.9, 8.1, 8.3, 12., 20.], cases=rows,
        pulse_admission='fixed_age_no_stationarity_gate',
        protocol='Critical groups: 2000 tau_c of full feedback after frozen preparation, then pulse. '
        'Observe >=500 tau_c after the longest pulse. tau_m=12/20 all1 retain original timing; '
        'reuse their eight completed controls. No drift/velocity veto on submission. '
        'Analyze signed same-seed pulse-control differences; retain finite-age interpretation.')
    (root / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (root / 'simulation_cases.txt').write_text(''.join(
        '/home/helu/mass_hd/cases/' + CAMPAIGN + '/' + r['case'] + '/run.dat\n'
        for r in rows if 'reuse_from' not in r))
    (root / 'README.md').write_text(
        '# Fixed-age pmem pulses\n\n'
        'Complete 160 previously held pulses: 128 near-critical all0/all1 pulses plus '
        '32 all1 pulses at tau_m/tau_c=12/20. Four seeds and all four pulse durations retained.\n\n'
        'Near-critical feedback waiting is 2000 tau_c after the original frozen preparation; '
        'recovery is at least 500 tau_c after the longest pulse. This is a prescribed age, '
        'not a guarantee of stationarity. No drift or velocity gate controls submission.\n\n'
        'Run 192 new simulations (32 controls +160 pulses), reuse eight existing 12/20 all1 '
        'controls with identical inputs. simulation_cases.txt excludes the reused controls. '
        'prepare_cw_pulse_reuse.py validates and reduces those raw archives without copying large videos.\n\n'
        'Summary uses cw_pmem_pulse_analysis.py --allow-drift. All paired curves remain available. '
        'A decaying paired response can be fitted despite baseline drift, retaining that diagnostic. '
        'A persistent late difference receives no forced zero-offset exponential fit.\n')
    print(json.dumps({'campaign': CAMPAIGN, 'new_simulations': 192, 'pulses': 160,
                      'new_controls': 32, 'reused_controls': 8,
                      'critical_nsteps': next(r['nsteps'] for r in rows if r['near_critical'])}))


if __name__ == '__main__':
    main()
