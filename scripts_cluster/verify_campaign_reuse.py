#!/usr/bin/env python3
"""Verify unchanged simulator sources and a previously verified SIF before reuse.

This is a provenance check, not a simulation or analysis. Run in the cluster
checkout after syncing the campaign commit. Derived evidence belongs on scratch.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def verify(manifest_path, previous_path, sif, expected_commit, out):
    repo = Path(__file__).resolve().parents[1]
    manifest = json.loads(manifest_path.read_text())
    previous = json.loads(previous_path.read_text())
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    if commit != expected_commit:
        raise ValueError('Local/cluster commit mismatch')
    source = manifest['simulator_source_sha256']
    if source != previous['source_sha256']:
        raise ValueError('Simulator sources differ from the verified image provenance')
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', 'HEAD',
                                     'src', 'Makefile'], cwd=repo, text=True).splitlines()
    if set(names) != set(source):
        raise ValueError('Tracked simulator source inventory changed')
    for name, digest in source.items():
        committed = subprocess.check_output(['git', 'show', f'HEAD:{name}'], cwd=repo)
        if hashlib.sha256(committed).hexdigest() != digest or sha(repo/name) != digest:
            raise ValueError(f'Cluster committed/worktree source mismatch: {name}')
    image_hash = sha(sif)
    # Accept only a byte-identical SIF, preserving the earlier embedded-source
    # and executable verification, despite subsequent case/analysis-only commits.
    if image_hash != previous['sif_sha256']:
        raise ValueError('SIF changed since previous source/binary verification')
    expected = set()
    for item in manifest['cases']:
        path = manifest_path.parent/item['case']/'run.dat'
        if sha(path) != item['run_dat_sha256']:
            raise ValueError(f'Runcard hash mismatch: {path}')
        expected.add(path.resolve())
    if expected != {p.resolve() for p in manifest_path.parent.rglob('*.dat')}:
        raise ValueError('Unexpected/missing .dat files in submission directory')
    result = {'campaign': manifest['campaign'], 'local_commit': expected_commit,
        'cluster_commit': commit, 'source_base_commit': manifest['source_base_commit'],
        'source_sha256': source, 'verified_source_files': len(source),
        'sif_path': str(sif), 'sif_resolved': str(sif.resolve()), 'sif_sha256': image_hash,
        'previous_provenance': str(previous_path), 'previous_provenance_sha256': sha(previous_path),
        'binary_sha256': previous['binary_sha256'],
        'case_manifest_sha256': sha(manifest_path), 'verified_cases': len(expected),
        'note': 'Reused byte-identical previously verified SIF; all committed and cluster '
                'worktree simulator sources match. Only campaign cases/analysis/workflow changed.'}
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise FileExistsError(f'Provenance already exists: {out}')
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('previous', type=Path)
    p.add_argument('sif', type=Path)
    p.add_argument('expected_commit')
    p.add_argument('out', type=Path)
    a = p.parse_args()
    verify(a.manifest, a.previous, a.sif, a.expected_commit, a.out)
