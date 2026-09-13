#!/usr/bin/env python3
"""Verify a newly built SIF, committed sources and campaign inputs before Slurm."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verify(manifest_path, sif, expected_commit, expected_binary):
    repo = Path(__file__).resolve().parents[1]
    manifest = json.loads(manifest_path.read_text())
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    if commit != expected_commit:
        raise ValueError('Local/cluster commit mismatch')
    source = manifest['simulator_source_sha256']
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    if set(names) != set(source):
        raise ValueError('Tracked simulator source inventory differs')
    runtime = ['scripts_cluster/submit_array.sh', 'scripts_cluster/submit_case.sh',
               'scripts_cluster/submit_analysis.sh', 'scripts_cluster/verify_campaign_image.py',
               'scripts_cluster/launch_cw_pmem_pulses.py',
               'plot/python/confluent_wet/cw_pmem_pulse_analysis.py']
    runtime_hashes = {}
    for name in names + runtime:
        committed = subprocess.check_output(['git', 'show', f'HEAD:{name}'], cwd=repo)
        digest = hashlib.sha256(committed).hexdigest()
        if sha(repo / name) != digest or (name in source and source[name] != digest):
            raise ValueError(f'Committed/worktree/manifest mismatch: {name}')
        if name in runtime:
            runtime_hashes[name] = digest
    inside = ['apptainer', 'exec', str(sif)]
    build = json.loads(subprocess.check_output(inside + ['cat', '/opt/mass_hd/build_info.json'], text=True))
    if build['commit'] != commit or build['source_sha256'] != source:
        raise ValueError('SIF build provenance does not match checkout and manifest')
    paths = ['/opt/mass_hd/' + name for name in names] + ['/opt/mass_hd/mass']
    actual = subprocess.check_output(inside + ['sha256sum'] + paths, text=True)
    hashes = {line.split(maxsplit=1)[1].strip(): line.split()[0] for line in actual.splitlines()}
    for name, digest in source.items():
        if hashes['/opt/mass_hd/' + name] != digest:
            raise ValueError(f'Actual SIF source mismatch: {name}')
    binary = hashes['/opt/mass_hd/mass']
    if binary != expected_binary:
        raise ValueError('SIF executable differs from the locally verified Docker binary')
    expected = set()
    for item in manifest['cases']:
        path = manifest_path.parent / item['case'] / 'run.dat'
        if sha(path) != item['run_dat_sha256']:
            raise ValueError(f'Runcard hash mismatch: {path}')
        relative = path.resolve().relative_to(repo)
        committed = subprocess.check_output(['git', 'show', f'HEAD:{relative}'], cwd=repo)
        if hashlib.sha256(committed).hexdigest() != item['run_dat_sha256']:
            raise ValueError(f'Runcard not committed as expected: {relative}')
        expected.add(path.resolve())
    if len(expected) != len(manifest['cases']) or expected != {p.resolve() for p in manifest_path.parent.rglob('*.dat')}:
        raise ValueError('Unexpected, missing or duplicate campaign inputs')
    return {'campaign': manifest['campaign'], 'local_commit': expected_commit,
            'cluster_commit': commit, 'image_build_commit': build['commit'],
            'source_sha256': source, 'runtime_sha256': runtime_hashes,
            'verified_source_files': len(source), 'verified_cases': len(expected),
            'sif_path': str(sif), 'sif_sha256': sha(sif), 'binary_sha256': binary,
            'case_manifest_sha256': sha(manifest_path)}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('sif', type=Path)
    p.add_argument('expected_commit')
    p.add_argument('expected_binary')
    p.add_argument('out', type=Path)
    a = p.parse_args()
    if a.out.exists():
        raise FileExistsError(a.out)
    result = verify(a.manifest, a.sif, a.expected_commit, a.expected_binary)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if not k.endswith('_sha256') or isinstance(v, str)}))
