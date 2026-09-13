#!/usr/bin/env python3
"""Generate paired controls and finite-duration pmem=0 pulses at mc=.2287."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

TC = 674.3290333006435
TM = (1., 2.5, 7.7, 7.9, 8.1, 8.3, 12., 20.)
TP = (.3, 1., 3., 10.)
CAMPAIGN = '20260913/cw_pmem_pulse'


def token(value):
    return f'{value:g}'.replace('.', 'p')


def generate(out):
    repo = Path(__file__).resolve().parents[1]
    baseline = repo/'cases/20260909/cw_mc_tm_phase/L256/mc0p21_tm20_chi0_rep1/run.dat'
    base = {}
    for line in baseline.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            k, v = line.split('=', 1)
            base[k.strip()] = v.strip()
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {n: hashlib.sha256(subprocess.check_output(['git', 'show', f'HEAD:{n}'], cwd=repo)).hexdigest()
              for n in names}
    rows = []
    for tm in TM:
        near = 7.7 <= tm <= 8.3
        prep = round(max(100., 10*tm)*TC)
        wait = round((500. if near else 100.)*TC)
        after = round((500. if near else 200.)*TC)
        start = prep+wait
        nsteps = start+round(max(TP)*TC)+after
        for init in (('patch',) if tm < 3 else ('0', '1')):
            for rep in range(1, 5):
                stem = f'L256/tm{token(tm)}_init{init}_rep{rep}'
                control = f'controls/{stem}'
                for tp in (0.,)+TP:
                    kind = 'control' if tp == 0 else 'pulse'
                    case = control if tp == 0 else f'pulses/tp{token(tp)}/{stem}'
                    p = dict(base)
                    chi0 = .5 if init == 'patch' else float(init)
                    m0 = .193 if init == 'patch' else .3301 if init == '0' else .0559
                    video_start = start-round(20*TC)
                    p.update({'mc': .2287, 'pmem': .016838, 'pmem-pulse-value': 0.,
                        'pmem-pulse-start': start, 'pmem-pulse-steps': round(tp*TC),
                        'nsteps': nsteps, 'ninfo': round(100*TC), 'nstart': 0,
                        'tau-m': tm*TC, 'tau-chi': 202.3, 'chi-freeze-steps': prep,
                        'mem-freeze-steps': 0, 'seed': 190900+rep, 'chi-seed': 200900+rep,
                        'chi-config': 'binary-noise' if init == 'patch' else 'uniform',
                        'chi0': chi0, 'm0': m0, 'chi-noise': 0.,
                        'chi-length': 7.02 if init == 'patch' else 0.,
                        'chi-lo': 0., 'chi-hi': 1., 'm-lo': .3301, 'm-hi': .0559,
                        'nresponse': 13, 'nvideo': round(TC), 'video-stride': 1,
                        'video-start': video_start, 'nvideo-dense': round(.05*TC),
                        'video-dense-start': video_start, 'video-dense-end': start+round(110*TC),
                        'frame-light': 1})
                    header = (f'# {CAMPAIGN}; paired temporary pmem reduction, no mechanical forcing.\n'
                        f'# tm/tc={tm:g}; tp/tc={tp:g}; tau_c={TC:.13f} steps.\n'
                        '# Pulse threshold 0; baseline .016838. Endpoints are a half-open step interval.\n'
                        f'# Preparation {prep}; full-feedback wait {wait}; pulse start {start}.\n'
                        f'# All five paired arms end at {nsteps}; recovery covers longest pulse plus {after} steps.\n'
                        '# Binary patches describe t=0, with matched m in each phase; preparation transports them.\n'
                        '# Controls are inspected before pulse submission; near-critical waiting is provisional.\n')
                    text = header+''.join(f'{k:<24} = {v:.15g}\n' if isinstance(v, float)
                                            else f'{k:<24} = {v}\n' for k, v in p.items())
                    path = out/case/'run.dat'
                    if path.exists() and path.read_text() != text:
                        raise FileExistsError(f'Refusing to replace differing runcard: {path}')
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(text)
                    rows.append({'case': case, 'kind': kind, 'paired_control': control,
                        'tm_over_tc': tm, 'tp_over_tc': tp, 'initialization': init,
                        'replicate': rep, 'seed': p['seed'], 'chi_seed': p['chi-seed'],
                        'L': 256, 'mc': .2287, 'preparation_steps': prep,
                        'feedback_wait_steps': wait, 'pulse_start_steps': start,
                        'pulse_duration_steps': round(tp*TC), 'recovery_steps_min': after,
                        'nsteps': nsteps, 'near_critical': near,
                        'run_dat_sha256': hashlib.sha256(text.encode()).hexdigest()})
    assert len(rows) == 280 and sum(r['kind'] == 'control' for r in rows) == 56
    manifest = {'campaign': CAMPAIGN, 'source_base_commit': commit,
        'simulator_source_sha256': hashes, 'baseline_runcard': str(baseline.relative_to(repo)),
        'fixed': {'mc': .2287, 'pmem': .016838, 'pmem_pulse_value': 0., 'tau_chi': 202.3,
                  'L': 256, 'r': .3, 'tau_c': TC},
        'tm_grid': TM, 'tp_grid': TP, 'replicates': 4,
        'protocol': 'Uniform pmem reduction in the memory source only; same-seed controls '
          'and pulse arms have identical prehistory. Controls retain the nominal pulse start. '
          'Inspect control stationarity before submitting pulses; extend only failed groups. '
          'All runs have enough time to compare 200/500 tau_c after the longest pulse.',
        'output': {'nresponse': 13, 'native_video_shape': [256, 256],
          'dense_video_dt_steps': round(.05*TC), 'sparse_video_dt_steps': round(TC),
          'dense_window_tc_relative_to_onset': [-20, 110],
          'note': 'Use video_meta.csv actual timestamps, including irregular event frames.'},
        'cases': rows}
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (out/'README.md').write_text(
        '# Temporary pmem=0 pulses at mc=.2287\n\n'
        'Eight memory times; binary patches at 1/2.5, all0/all1 at other times; four independent seeds. '
        '224 pulse runs plus 56 controls. Pulse durations .3/1/3/10 reference tau_c.\n\n'
        'Baseline pmem=.016838, pulse pmem=0. No density reset or direct flow forcing. '
        'Controls and all pulse durations repeat the same seed and history.\n\n'
        'Feedback waiting: 100 tau_c away from the transition; 500 tau_c at 7.7/7.9/8.1/8.3. '
        'Preceded by the original max(100 tau_c,10 tau_m) reaction-frozen preparation. '
        'Post-pulse observation: at least 200/500 tau_c respectively. '
        'Short near-critical preparation is a first control check, not a stationarity guarantee.\n\n'
        'response.csv samples full-grid moments and source occupancies every 13 steps and at events. '
        'Video stores all 256x256 sites from 20 tau_c before onset; dense through onset+110 tau_c, '
        'then sparse. Actual timestamps must determine movie speed. '
        'Late exponential relaxation is fitted only to resolved return-to-branch responses.\n\n'
        'Submit controls with cw_pmem_pulse_analysis.py as the per-case reducer. '
        'Then run scripts_cluster/launch_cw_pmem_pulses.py as an afterany analysis job. '
        'It requires all four seeds of each tau_m/initialization group to pass the final '
        '50/200 tau_c baseline and control-tail drift checks before submitting the 16 paired pulses. '
        'Failed groups remain explicitly held in control_gate/gate.json for locally generated '
        'longer preparation cases; individual successful seeds are not selected in isolation.\n')
    print(json.dumps({'out': str(out), 'controls': 56, 'pulses': 224,
                      'nsteps_min': min(r['nsteps'] for r in rows),
                      'nsteps_max': max(r['nsteps'] for r in rows)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out', type=Path)
    generate(parser.parse_args().out)
