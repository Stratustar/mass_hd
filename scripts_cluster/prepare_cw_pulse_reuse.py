#!/usr/bin/env python3
"""Validate/reduce completed controls reused by a new pulse campaign. Slurm only."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


def main(manifest_path, poster=False):
    repo = Path(__file__).resolve().parents[1]
    manifest = json.loads(manifest_path.read_text())
    scratch = Path('/scratch/helu/mass_hd')
    raw = scratch / 'cases' / manifest['campaign']
    results = scratch / 'results/cases' / manifest['campaign']
    if poster:
        sys.path.insert(0, str(repo/'plot/python/confluent_wet'))
        import cw_poster_pulse as analysis
        analysis.load_manifest(manifest_path)
        original_manifest = repo/'cases'/manifest['source_campaign']/'manifest.json'
        if analysis.sha(original_manifest) != manifest['source_manifest_sha256']:
            raise ValueError('Source control campaign manifest changed')
    else:
        spec = importlib.util.spec_from_file_location('pulse_analysis',
            repo / 'plot/python/confluent_wet/cw_pmem_pulse_analysis.py')
        analysis = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(analysis)
    reused = []
    for item in manifest['cases']:
        if 'reuse_from' not in item:
            continue
        assert item['kind'] == 'control'
        original_card = repo / 'cases' / item['reuse_from'] / 'run.dat'
        new_card = manifest_path.parent / item['case'] / 'run.dat'
        if (original_card.read_bytes() != new_card.read_bytes() or
                hashlib.sha256(new_card.read_bytes()).hexdigest() != item['run_dat_sha256']):
            raise ValueError('Reused control does not have identical inputs')
        source = scratch / 'cases' / item['reuse_from']
        destination = raw / item['case']
        if poster:
            # Re-reduce the completed raw control with current strict auditing;
            # do not relabel an old analysis provenance record as a new one.
            analysis.reduce_case(source, results/item['case'], manifest_path)
            row = json.loads((results/item['case']/'pulse_case.json').read_text())
            analysis.load_case(results/item['case'], item, manifest, manifest_path)
        else:
            row = analysis.reduce_case(source, results / item['case'])
        for key in ('nsteps', 'preparation_steps', 'pulse_start_steps', 'pulse_duration_steps',
                    'initialization', 'replicate', 'seed'):
            if row[key] != item[key]:
                raise ValueError('Reused control differs in ' + key)
        if destination.exists() or destination.is_symlink():
            if not destination.is_symlink() or destination.resolve() != source.resolve():
                raise FileExistsError(destination)
        else:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.symlink_to(source, target_is_directory=True)
        reused.append({'case': item['case'], 'source': str(source),
                       'input_sha256': item['run_dat_sha256'],
                       'response_sha256': hashlib.sha256((source / 'response.csv').read_bytes()).hexdigest()})
    if poster and len(reused) != manifest['reused_controls']:
        raise ValueError('Wrong reused-control count')
    results.mkdir(parents=True, exist_ok=True)
    (results / 'reused_controls.json').write_text(json.dumps(reused, indent=2) + '\n')
    print(json.dumps({'validated_reused_controls': len(reused)}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('--poster', action='store_true')
    args = p.parse_args()
    main(args.manifest, poster=args.poster)
