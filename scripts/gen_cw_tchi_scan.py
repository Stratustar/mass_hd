#!/usr/bin/env python3
"""Phenotype-clock scan at mc=0.2287, tau_m=0.3 tau_c: vary tau_chi, uniform 0/1 starts.

Each input is copied from the native-resolution tm0p3 uniform runcard of
20260915/cw_mc02287_four_init_fullres (seed rep1). Only tau-chi and nsteps change:
the observation is shortened from 2000 to 1000 tau_c, all physics, seeds,
preparation and native 256x256 video output are retained.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

TC = 674.3290333006435
TM = 0.3
TCHI = (0.3, 0.6, 1.2, 2.5, 5.0, 10.0, 20.0)
STARTS = ('0', '1')
PREP = 67433                         # round(100 tau_c), as in the source campaign
OBS = math.ceil(1000*TC)             # >= 1000 tau_c so the 500 tau_c tail fits twice
SOURCE = 'cases/20260915/cw_mc02287_four_init_fullres'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def parameters(text):
    return dict(tuple(x.strip() for x in line.split('=', 1))
                for line in text.splitlines()
                if line.strip() and not line.lstrip().startswith('#'))


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
    source_manifest = repo/SOURCE/'manifest.json'
    src_rows = {(r['initialization'], r['replicate']): r
                for r in json.loads(source_manifest.read_text())['cases']
                if r['tm_over_tc'] == TM}
    rows, pending = [], []
    for start in STARTS:
        src = src_rows[(start, 1)]
        source = repo/SOURCE/src['case']/'run.dat'
        old = source.read_bytes()
        if digest(old) != src['run_dat_sha256']:
            raise ValueError(f'Source runcard changed: {source}')
        old_text = old.decode()
        before = parameters(old_text)
        assert before['LX'] == before['LY'] == '256' and before['video-stride'] == '1'
        assert before['nvideo'] == '337' and before['chi-config'] == 'uniform'
        assert int(before['chi-freeze-steps']) == PREP and before['mc'] == '0.2287'
        assert math.isclose(float(before['tau-m'])/TC, TM, rel_tol=1e-9)
        body = old_text.split('\nmodel', 1)
        assert len(body) == 2
        for tc in TCHI:
            tau_chi = round(tc*TC, 10)
            text = 'model' + body[1]
            text = replace(text, 'tau-chi', f'{tau_chi:.10g}')
            text = replace(text, 'nsteps', PREP+OBS)
            header = (
                f'# {campaign}: phenotype-clock scan at mc=0.2287, tau_m/tau_c={TM:g}.\n'
                f'# tau_chi/tau_c={tc:g} ({tau_chi:.10g} steps); tau_c={TC!r} steps.\n'
                f'# Uniform start chi={start}; seed rep1. Copied from {src["case"]} of {SOURCE};\n'
                '# only tau-chi and nsteps differ (observation 2000 -> 1000 tau_c).\n'
                f'# Prepare {PREP} steps with chi REACTION off (transport on); then observe {OBS} steps.\n'
                '# Native 256x256 u/P/m/chi video every 337 steps; exact full-field means at the same cadence.\n')
            text = header + text
            after = parameters(text)
            changed = {k for k in before.keys() | after.keys() if before.get(k) != after.get(k)}
            assert changed == ({'nsteps'} if after['tau-chi'] == before['tau-chi']
                               else {'nsteps', 'tau-chi'}), changed
            case = f'L256/mc0p2287_tc{token(tc)}_tm0p3_chi{start}_rep1'
            rows.append({'case': case, 'L': 256, 'mc': 0.2287, 'tm_over_tc': TM,
                         'tau_chi_over_tc': tc, 'tau_chi_steps': tau_chi,
                         'initialization': start, 'chi0': float(start), 'replicate': 1,
                         'seed': int(after['seed']), 'chi_seed': int(after['chi-seed']),
                         'preparation_steps': PREP, 'observation_steps': OBS,
                         'nsteps': PREP+OBS, 'run_dat_sha256': digest(text.encode()),
                         'source_case': str(source.relative_to(repo)),
                         'source_run_dat_sha256': digest(old)})
            pending.append((out/case/'run.dat', text))
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {name: digest(subprocess.check_output(['git', 'show', f'HEAD:{name}'], cwd=repo))
              for name in names}
    frames = (PREP+OBS)//337+1
    manifest = {
        'campaign': campaign, 'tau_c_reference': TC, 'source_campaign': SOURCE,
        'source_manifest_sha256': digest(source_manifest.read_bytes()),
        'source_base_commit': commit, 'simulator_source_sha256': hashes,
        'fixed': {'L': 256, 'mc': 0.2287, 'pmem': 0.016838, 'r': 0.3,
                  'tau_m_over_tc': TM, 'tau_m_steps': float(parameters(pending[0][1])['tau-m'])},
        'tchi_grid': list(TCHI), 'initializations': list(STARTS), 'replicates': [1],
        'changed_parameters': {'tau-chi': 'tau_chi/tau_c * tau_c, rounded to 1e-10 steps '
                               '(the source used 202.3 for 0.3 tau_c; here 202.2987099902)',
                               'nsteps': f'{PREP}+{OBS} (observation 1000 tau_c, was 2000)'},
        'protocol': 'Preparation 100 tau_c with chi reaction frozen, flow/memory/advection/'
                    'diffusion running; observe 1000 tau_c; terminal window 500 tau_c.',
        'statistics': {'observable': 'exact full-domain mean chi', 'tail_tc': 500},
        'output': {'native_shape': [256, 256], 'video_stride': 1, 'nvideo': 337,
                   'ninfo': 67400, 'frame_light': 1, 'stored_fields': ['u', 'P', 'm', 'chi'],
                   'frames_per_run': frames, 'raw_video_bytes': len(rows)*frames*4*256*256},
        'executable_note': 'No executable change; reuse the verified simulator SIF.',
        'cases': rows}
    pending += [(out/'manifest.json', json.dumps(manifest, indent=2)+'\n'),
                (out/'README.md',
                 '# tau_chi scan at mc=0.2287, tau_m=0.3 tau_c\n\n'
                 f'{len(rows)} runs: tau_chi/tau_c = {", ".join(f"{x:g}" for x in TCHI)}; '
                 'uniform chi=0 (m0=.3301) and chi=1 (m0=.0559); one seed (190901/200901).\n\n'
                 f'Copied from the tm0p3 rep1 uniform cards of {SOURCE}; only tau-chi and nsteps '
                 'change. L256, mc=.2287, pmem=.016838, r=.3, tau_m=202.2987 steps retained.\n\n'
                 'Preparation 100 tau_c (chi reaction off), observation 1000 tau_c, '
                 f'{PREP+OBS} steps. Native 256x256 uint8 u/P/m/chi every 337 steps '
                 f'(~0.5 tau_c), {frames} frames/run.\n\n'
                 'Per-run analysis and summary: plot/python/confluent_wet/cw_tchi_scan.py.\n')]
    for path, body in pending:
        if path.exists() and path.read_text() != body:
            raise FileExistsError(f'Refusing to replace differing file: {path}')
    for path, body in pending:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body)
    assert {p.resolve() for p in out.rglob('*.dat')} == {
        (out/r['case']/'run.dat').resolve() for r in rows}
    print(json.dumps({'campaign': campaign, 'cases': len(rows), 'tchi_grid': TCHI,
                      'nsteps': PREP+OBS, 'output': manifest['output']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    generate(p.parse_args().out)
