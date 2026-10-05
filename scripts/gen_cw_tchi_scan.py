#!/usr/bin/env python3
"""Phenotype-clock scan at mc=0.2287 and one memory time: vary tau_chi, uniform 0/1 starts.

Each input is copied from the native-resolution uniform runcard (seed rep1) of
20260915/cw_mc02287_four_init_fullres at the requested tau_m grid coordinate. Only
tau-chi and, optionally, nsteps change (--obs-tc shortens the 2000 tau_c observation);
all physics, seeds, preparation and native 256x256 video output are retained.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

TC = 674.3290333006435
TCHI = (0.3, 0.6, 1.2, 2.5, 5.0, 10.0, 20.0)
STARTS = ('0', '1')
PREP = 67433                         # round(100 tau_c), as in the source campaign
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


def generate(out, TM, obs_tc):
    # ceil keeps the observation >= obs_tc, so the reducer's 500 tau_c tail fits twice.
    OBS = math.ceil(obs_tc*TC) if obs_tc else None
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
        obs = OBS or int(before['nsteps'])-PREP
        tm_tok = re.search(r'_tm([0-9p]+)_chi', src['case']).group(1)
        body = old_text.split('\nmodel', 1)
        assert len(body) == 2
        for tc in TCHI:
            tau_chi = round(tc*TC, 10)
            text = 'model' + body[1]
            text = replace(text, 'tau-chi', f'{tau_chi:.10g}')
            text = replace(text, 'nsteps', PREP+obs)
            header = (
                f'# {campaign}: phenotype-clock scan at mc=0.2287, tau_m/tau_c={TM!r}.\n'
                f'# tau_chi/tau_c={tc:g} ({tau_chi:.10g} steps); tau_c={TC!r} steps.\n'
                f'# Uniform start chi={start}; seed rep1. Copied from {src["case"]} of {SOURCE};\n'
                f'# only tau-chi{" and nsteps differ" if OBS else " differs"} (observation {obs/TC:.1f} tau_c).\n'
                f'# Prepare {PREP} steps with chi REACTION off (transport on); then observe {obs} steps.\n'
                '# Native 256x256 u/P/m/chi video every 337 steps; exact full-field means at the same cadence.\n')
            text = header + text
            after = parameters(text)
            changed = {k for k in before.keys() | after.keys() if before.get(k) != after.get(k)}
            assert changed - {'nsteps', 'tau-chi'} == set() and ('nsteps' in changed) == bool(OBS), changed
            case = f'L256/mc0p2287_tc{token(tc)}_tm{tm_tok}_chi{start}_rep1'
            rows.append({'case': case, 'L': 256, 'mc': 0.2287, 'tm_over_tc': TM,
                         'tau_chi_over_tc': tc, 'tau_chi_steps': tau_chi,
                         'initialization': start, 'chi0': float(start), 'replicate': 1,
                         'seed': int(after['seed']), 'chi_seed': int(after['chi-seed']),
                         'preparation_steps': PREP, 'observation_steps': obs,
                         'nsteps': PREP+obs, 'run_dat_sha256': digest(text.encode()),
                         'source_case': str(source.relative_to(repo)),
                         'source_run_dat_sha256': digest(old)})
            pending.append((out/case/'run.dat', text))
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {name: digest(subprocess.check_output(['git', 'show', f'HEAD:{name}'], cwd=repo))
              for name in names}
    obs = rows[0]['observation_steps']
    frames = (PREP+obs)//337+1
    manifest = {
        'campaign': campaign, 'tau_c_reference': TC, 'source_campaign': SOURCE,
        'source_manifest_sha256': digest(source_manifest.read_bytes()),
        'source_base_commit': commit, 'simulator_source_sha256': hashes,
        'fixed': {'L': 256, 'mc': 0.2287, 'pmem': 0.016838, 'r': 0.3,
                  'tau_m_over_tc': TM, 'tau_m_steps': float(parameters(pending[0][1])['tau-m'])},
        'tchi_grid': list(TCHI), 'initializations': list(STARTS), 'replicates': [1],
        'changed_parameters': {'tau-chi': 'tau_chi/tau_c * tau_c, rounded to 1e-10 steps '
                               '(the source used 202.3 for 0.3 tau_c; here 202.2987099902)',
                               'nsteps': (f'{PREP}+{obs} (observation {obs/TC:.1f} tau_c, was 2000)'
                                          if OBS else 'unchanged (observation 2000 tau_c)')},
        'protocol': 'Preparation 100 tau_c with chi reaction frozen, flow/memory/advection/'
                    f'diffusion running; observe {obs/TC:.1f} tau_c; terminal window 500 tau_c.',
        'statistics': {'observable': 'exact full-domain mean chi', 'tail_tc': 500},
        'output': {'native_shape': [256, 256], 'video_stride': 1, 'nvideo': 337,
                   'ninfo': 67400, 'frame_light': 1, 'stored_fields': ['u', 'P', 'm', 'chi'],
                   'frames_per_run': frames, 'raw_video_bytes': len(rows)*frames*4*256*256},
        'executable_note': 'No executable change; reuse the verified simulator SIF.',
        'cases': rows}
    pending += [(out/'manifest.json', json.dumps(manifest, indent=2)+'\n'),
                (out/'README.md',
                 f'# tau_chi scan at mc=0.2287, tau_m={TM!r} tau_c\n\n'
                 f'{len(rows)} runs: tau_chi/tau_c = {", ".join(f"{x:g}" for x in TCHI)}; '
                 'uniform chi=0 (m0=.3301) and chi=1 (m0=.0559); one seed (190901/200901).\n\n'
                 f'Copied from the tm{tm_tok} rep1 uniform cards of {SOURCE}; only tau-chi'
                 f'{" and nsteps" if OBS else ""} change. L256, mc=.2287, pmem=.016838, r=.3, '
                 f'tau_m={manifest["fixed"]["tau_m_steps"]:g} steps retained.\n\n'
                 f'Preparation 100 tau_c (chi reaction off), observation {obs/TC:.1f} tau_c, '
                 f'{PREP+obs} steps. Native 256x256 uint8 u/P/m/chi every 337 steps '
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
                      'nsteps': PREP+obs, 'output': manifest['output']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    p.add_argument('--tm', type=float, default=0.3, help='tau_m/tau_c, an exact source grid coordinate')
    p.add_argument('--obs-tc', type=float, default=None,
                   help='observation in tau_c (default: keep the source 2000 tau_c)')
    a = p.parse_args()
    generate(a.out, a.tm, a.obs_tc)
