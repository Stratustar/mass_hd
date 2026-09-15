#!/usr/bin/env python3
"""Copy five original four-start scan groups, changing only video-stride to 1."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import subprocess

TM = (0.3, 6.821052632, 7.752631579, 8.684210526, 18.0)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def parameters(text):
    return dict(tuple(x.strip() for x in line.split('=', 1))
                for line in text.splitlines()
                if line.strip() and not line.lstrip().startswith('#'))


def generate(out):
    repo = Path(__file__).resolve().parents[1]
    out = out.resolve()
    campaign = out.relative_to(repo / 'cases').as_posix()
    original = repo / 'cases/20260911/cw_mc02287_four_init'
    source_manifest = original / 'manifest.json'
    manifest = json.loads(source_manifest.read_text())
    rows = [copy.deepcopy(r) for r in manifest['cases'] if r['tm_over_tc'] in TM]
    expected = {(t, s, r) for t in TM for s in ('0', '1', 'stripe01', 'noise')
                for r in (1, 2)}
    assert len(rows) == 40
    assert {(r['tm_over_tc'], r['initialization'], r['replicate']) for r in rows} == expected
    pending = []
    for row in rows:
        source = original / row['case'] / 'run.dat'
        old = source.read_bytes()
        if digest(old) != row['run_dat_sha256']:
            raise ValueError(f'Original runcard changed: {source}')
        old_text = old.decode()
        before = parameters(old_text)
        assert before['video-stride'] == '8'
        assert before['LX'] == before['LY'] == '256' and before['nvideo'] == '337'
        body, count = re.subn(r'(?m)^video-stride\s*=\s*8\s*$',
                             'video-stride         = 1', old_text)
        assert count == 1
        body = (f'# Full-resolution rerun: {campaign}; original input below.\n'
                '# Only video-stride changes (8 -> 1); all physics, seeds and clocks retained.\n'
                + body)
        after = parameters(body)
        assert {k for k in before.keys() | after.keys()
                if before.get(k) != after.get(k)} == {'video-stride'}
        assert after['video-stride'] == '1'
        row.update(source_case=str(source.relative_to(repo)),
                   source_run_dat_sha256=digest(old), run_dat_sha256=digest(body.encode()))
        pending.append((out / row['case'] / 'run.dat', body))
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    hashes = {name: digest(subprocess.check_output(['git', 'show', f'HEAD:{name}'], cwd=repo))
              for name in names}
    manifest.update(campaign=campaign, source_base_commit=commit,
                    simulator_source_sha256=hashes, tm_grid=list(TM), cases=rows,
                    tm_grid_rule='Five exact existing grid coordinates selected for full-resolution videos',
                    original_campaign='20260911/cw_mc02287_four_init',
                    original_manifest_sha256=digest(source_manifest.read_bytes()),
                    changed_parameters={'video-stride': {'original': 8, 'new': 1}},
                    output={'native_shape': [256, 256], 'video_stride': 1, 'nvideo': 337,
                            'ninfo': 67400, 'frame_light': 1,
                            'stored_fields': ['u', 'P', 'm', 'chi'], 'quantization': 'uint8',
                            'total_video_frames': sum(r['nsteps']//337+1 for r in rows),
                            'raw_video_bytes': sum(r['nsteps']//337+1 for r in rows)*4*256*256},
                    executable_note='Use current verified simulator with pulse, response and windowed-video '
                                    'features disabled by their defaults. No executable changes in this campaign.')
    pending += [(out / 'manifest.json', json.dumps(manifest, indent=2)+'\n'),
                (out / 'README.md',
                 '# Full-resolution videos: five four-initialization groups\n\n'
                 '40 runs: tau_m/tau_c = 0.3, 6.821052632, 7.752631579, 8.684210526, 18; '
                 'all0, all1, stripe01 and Gaussian noise; original paired seeds 190901/190902.\n\n'
                 'Each input is copied from the original 20260911 campaign with only '
                 'video-stride changed from 8 to 1. L256, mc=0.2287, pmem=0.016838, '
                 'tau_chi=202.3 steps and all other dynamics are retained. No pressure pulse.\n\n'
                 'Preparation 100–180 tau_c; feedback observation 2000 tau_c; '
                 'total 1,416,091–1,470,037 steps. The original t=0 mixed initializations '
                 'continue to advect/diffuse during preparation.\n\n'
                 'Four native256x256 uint8 streams every337steps (~0.499756tau_c), '
                 '169400 total stored timestamps,44.407GB raw streams. Exact full-domain '
                 'scalar statistics retain their existing cadence. Original sparse JSON snapshots '
                 'and fixed stored value ranges are retained.\n\n'
                 'Per-run analysis: cw_four_init_scan.py. Summary: same script with '
                 '--summary --manifest. New output directories keep original results intact.\n')]
    for path, body in pending:
        if path.exists() and path.read_text() != body:
            raise FileExistsError(f'Refusing to replace differing file: {path}')
    for path, body in pending:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body)
    assert {p.resolve() for p in out.rglob('*.dat')} == {
        (out/r['case']/'run.dat').resolve() for r in rows}
    print(json.dumps({'campaign': campaign, 'cases': len(rows), 'tau_m_over_tau_c': TM,
                      'changed_parameters': manifest['changed_parameters'],
                      'output': manifest['output']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    generate(p.parse_args().out)
