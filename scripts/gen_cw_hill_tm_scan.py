#!/usr/bin/env python3
"""Memory-time scan of the HILL response law at mc=0.2287, beta=0: uniform active/passive starts.

Each input is copied from the rep1 uniform runcard of 20260911/cw_mc02287_four_init at the
same tau_m grid coordinate, so every point has a legacy step-law partner run. Changes:

  response law  switch-model=hill, s(P) logistic (wp), h(m) Hill (h-n), beta=0
  convention    chi=1 is ACTIVE under the hill law, so the legacy all-0 (active) start
                becomes chi0=1 and the legacy all-1 (passive) start chi0=0; m0 is unchanged
  diffusivity   Dm = Dchi = the legacy Dbio, Dbio = 0 (bit-identical to Dbio, checked by
                scripts/check_cw_hill.py)
  video         native 256x256 (video-stride 8 -> 1), as in 20260915 fullres

Everything else -- seeds, preparation, observation, clocks, pmem, mc, tau_chi -- is retained,
except that --tchi-over-tm R sets tau_chi = R tau_m point by point.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

TC = 674.3290333006435
TM = (0.3, 2.163157895, 4.957894737, 7.752631579, 8.684210526, 12.410526316, 18.0)
SOURCE = 'cases/20260911/cw_mc02287_four_init'
HILL = {'h-n': 16, 'wp': .002, 'beta': 0, 'g-n': 2, 'g-chic': .5, 'g-lock': 'active'}
DBIO = '0.00166667'
# legacy start label -> hill chi0 (the m0 of the source card is kept)
FLIP = {'0': '1', '1': '0'}


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


def generate(out, tchi_over_tm=None):
    """tchi_over_tm: None keeps the source tau_chi (202.3 steps = 0.3 tau_c) at every tau_m;
    a number R sets tau_chi = R tau_m point by point (the value tau-m carries in the card)."""
    repo = Path(__file__).resolve().parents[1]
    out = out.resolve()
    campaign = out.relative_to(repo/'cases').as_posix()
    source_manifest = repo/SOURCE/'manifest.json'
    src_rows = {(r['tm_over_tc'], r['initialization']): r
                for r in json.loads(source_manifest.read_text())['cases']
                if r['replicate'] == 1 and r['initialization'] in FLIP}
    variant = f'hn{token(HILL["h-n"])}_wp{token(HILL["wp"])}_b{token(HILL["beta"])}'
    rows, pending = [], []
    for tm in TM:
        for legacy_start, start in FLIP.items():
            src = src_rows[(tm, legacy_start)]
            source = repo/SOURCE/src['case']/'run.dat'
            old = source.read_bytes()
            if digest(old) != src['run_dat_sha256']:
                raise ValueError(f'Source runcard changed: {source}')
            old_text = old.decode()
            before = parameters(old_text)
            assert before['LX'] == before['LY'] == '256' and before['video-stride'] == '8'
            assert before['chi-config'] == 'uniform' and before['chi0'] == legacy_start
            assert before['switch-sign'] == '-1' and before['chi-width'] == '0'
            assert before['pmem-width'] == '0' and before['Dbio'] == DBIO
            assert before['mc'] == '0.2287' and abs(float(before['tau-m'])/TC/tm - 1) < 1e-9
            tm_tok = re.search(r'_tm([0-9p]+)_chi', src['case']).group(1)
            text = 'model' + old_text.split('\nmodel', 1)[1]
            for key, value in (('video-stride', 1), ('switch-sign', 1), ('chi0', start),
                               ('Dbio', 0)):
                text = replace(text, key, value)
            if tchi_over_tm is not None:
                tau_chi = f'{tchi_over_tm*float(before["tau-m"]):.13g}'
                text = replace(text, 'tau-chi', tau_chi)
            tchi_note = (f'tau_chi = {tchi_over_tm:g} tau_m = {tau_chi} steps'
                         if tchi_over_tm is not None else
                         f'tau_chi = {before["tau-chi"]} steps (source value)')
            text += ''.join(f'{k:<21}= {v}\n' for k, v in
                            [('switch-model', 'hill'), ('Dm', DBIO), ('Dchi', DBIO)]
                            + list(HILL.items()))
            header = (
                f'# {campaign}: hill response law, beta=0, at mc=0.2287, tau_m/tau_c={tm!r}.\n'
                f'# h(m)=m^n/(mc^n+m^n), n={HILL["h-n"]}; s(P)=1/(1+exp(-(P-pmem)/wp)), '
                f'wp={HILL["wp"]}; no phenotype lock.\n'
                '# chi=1 is ACTIVE here (zeta[z0+(1-z0)chi]); the legacy start '
                f'chi={legacy_start} is chi={start}, same m0.\n'
                f'# Copied from {src["case"]} of {SOURCE}; Dm=Dchi=old Dbio; '
                'native 256x256 video. Seeds, clocks, preparation retained.\n'
                + (f'# {tchi_note}.\n' if tchi_over_tm is not None else ''))
            text = header + text
            after = parameters(text)
            changed = {k for k in before.keys() | after.keys() if before.get(k) != after.get(k)}
            assert changed == {'video-stride', 'switch-sign', 'chi0', 'Dbio', 'switch-model',
                               'Dm', 'Dchi'} | set(HILL) | (
                               {'tau-chi'} if tchi_over_tm is not None else set()), changed
            tc_tok = f'_tchi{token(tchi_over_tm)}tm' if tchi_over_tm is not None else ''
            case = f'L256/mc0p2287_{variant}{tc_tok}_tm{tm_tok}_chi{start}_rep1'
            rows.append({'case': case, 'L': 256, 'mc': 0.2287, 'tm_over_tc': tm,
                         'initialization': start, 'chi0': float(start),
                         'legacy_initialization': legacy_start, 'replicate': 1,
                         'seed': int(after['seed']), 'chi_seed': int(after['chi-seed']),
                         'm0': float(after['m0']),
                         'preparation_steps': int(after['chi-freeze-steps']),
                         'nsteps': int(after['nsteps']), 'tau_chi_steps': float(after['tau-chi']),
                         'run_dat_sha256': digest(text.encode()),
                         'source_case': str(source.relative_to(repo)),
                         'source_run_dat_sha256': digest(old)})
            pending.append((out/case/'run.dat', text))
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    # The cards and the hill law are committed TOGETHER, so the simulator inventory is read
    # from the INDEX: the staged hill law, and the committed version of anything with
    # unstaged worktree edits (which are not part of the build). The hill law itself must
    # be staged exactly as it is on disk.
    hill_src = ['src/models/confluent-wet.cpp', 'src/models/confluent-wet.hpp']
    if subprocess.run(['git', 'diff', '--quiet', '--'] + hill_src, cwd=repo).returncode:
        raise ValueError('Stage the hill-law sources before generating the campaign')
    names = subprocess.check_output(['git', 'ls-files', 'src', 'Makefile'], cwd=repo,
                                    text=True).splitlines()
    hashes = {name: digest(subprocess.check_output(['git', 'show', f':{name}'], cwd=repo))
              for name in names}
    manifest = {
        'campaign': campaign, 'tau_c_reference': TC, 'source_campaign': SOURCE,
        'source_manifest_sha256': digest(source_manifest.read_bytes()),
        'generator_base_commit': commit, 'source_base_commit': commit,
        'simulator_source_sha256': hashes,
        'fixed': {'L': 256, 'mc': 0.2287, 'pmem': 0.016838, 'zeta0_frac': 0.3,
                  **({'tau_chi_over_tm': tchi_over_tm} if tchi_over_tm is not None
                     else {'tau_chi_steps': 202.3}),
                  'Dm': float(DBIO), 'Dchi': float(DBIO), **HILL},
        'tm_grid': list(TM), 'initializations': ['0', '1'], 'replicates': [1],
        'convention': 'switch-model=hill: chi=1 ACTIVE. chi0=1 is the legacy all-0 start.',
        'changed_parameters': ['video-stride 8->1', 'switch-sign -1->1', 'chi0 flipped',
                               'Dbio->0 with Dm=Dchi=old Dbio', 'switch-model=hill'] +
                              [f'{k}={v}' for k, v in HILL.items()],
        'protocol': 'Source protocol retained: preparation max(100 tau_c, ...) with chi reaction '
                    'frozen, observation 2000 tau_c, terminal window max(500, 20 tau_m) tau_c.',
        'output': {'native_shape': [256, 256], 'video_stride': 1, 'nvideo': 337,
                   'ninfo': 67400, 'frame_light': 1},
        'cases': rows}
    pending += [(out/'manifest.json', json.dumps(manifest, indent=2)+'\n'),
                (out/'README.md',
                 '# Hill response law, beta = 0: tau_m scan at mc = 0.2287\n\n'
                 f'{len(rows)} runs: tau_m/tau_c = {", ".join(f"{x:g}" for x in TM)}; '
                 + (f'tau_chi = {tchi_over_tm:g} tau_m at every point; ' if tchi_over_tm is not None
                    else 'tau_chi = 202.3 steps (0.3 tau_c); ') + 'uniform '
                 'chi=1 (active, m0=.3301) and chi=0 (passive, m0=.0559); one seed '
                 '(190901/200901).\n\n'
                 f'h-n={HILL["h-n"]} (tanh-equivalent width 2mc/n = {2*.2287/HILL["h-n"]:.4f}), '
                 f'wp={HILL["wp"]} (~0.1 sigma_P), beta=0, Dm=Dchi={DBIO}. Each card is the '
                 f'rep1 uniform card of {SOURCE} at the same tau_m with the response law '
                 'swapped and chi flipped (chi=1 ACTIVE), so every point has a legacy partner.\n\n'
                 'Per-run analysis and summary: plot/python/confluent_wet/cw_hill_tm_scan.py.\n')]
    for path, body in pending:
        if path.exists() and path.read_text() != body:
            raise FileExistsError(f'Refusing to replace differing file: {path}')
    for path, body in pending:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body)
    assert {p.resolve() for p in out.rglob('*.dat')} == {
        (out/r['case']/'run.dat').resolve() for r in rows}
    print(json.dumps({'campaign': campaign, 'cases': len(rows), 'tm_grid': TM,
                      'nsteps': sorted({r['nsteps'] for r in rows})}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    p.add_argument('--tchi-over-tm', type=float, default=None,
                   help='set tau_chi = R tau_m at every point (default: keep 202.3 steps)')
    a = p.parse_args()
    generate(a.out, a.tchi_over_tm)
