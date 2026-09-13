#!/usr/bin/env python3
"""After controls finish, screen whole four-seed groups and submit eligible pulses.

Run through submit_analysis.sh with afterany:<control-array>. Failed groups are
held explicitly for a longer preparation protocol; no seed is silently dropped.
All paths selected here refer to committed inputs. Derived files stay on scratch.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess

from verify_campaign_image import verify


def choose_groups(manifest, summary):
    controls = [r for r in manifest['cases'] if r['kind'] == 'control']
    groups = {}
    for row in controls:
        groups.setdefault((row['tm_over_tc'], row['initialization']), []).append(row)
    report, eligible = [], set()
    for (tm, init), rows in sorted(groups.items()):
        reasons = {}
        if {r['replicate'] for r in rows} != {1, 2, 3, 4} or len(rows) != 4:
            raise ValueError('Each initialization group must contain four independent controls')
        for item in rows:
            record = summary['cases'].get(item['case'])
            if record is None:
                reasons[item['case']] = ['missing_or_invalid_control']
            elif not record['control_usable']:
                reasons[item['case']] = record['baseline']['flags'] + record['control_post_flags']
        if not reasons:
            eligible.update(r['case'] for r in rows)
        report.append({'tm_over_tc': tm, 'initialization': init,
                       'eligible': not reasons, 'control_failures': reasons})
    selected = [r for r in manifest['cases'] if r['kind'] == 'pulse' and r['paired_control'] in eligible]
    return report, selected


def job_id(output):
    found = re.findall(r'^(?:Submitted batch job )?(\d+)(?:;\S+)?$', output, re.M)
    if len(found) != 1:
        raise ValueError('Could not identify exactly one submitted Slurm job: ' + output)
    return found[0]


def launch(manifest_path, results, sif_name, provenance_path):
    repo = Path(__file__).resolve().parents[1]
    manifest_path = manifest_path.resolve()
    results = results.resolve()
    results.relative_to(Path('/scratch/helu/mass_hd/results'))
    provenance = json.loads(provenance_path.read_text())
    # Recheck at dispatch time: another campaign may have updated the checkout.
    current = verify(manifest_path, Path(provenance['sif_path']),
                     provenance['local_commit'], provenance['binary_sha256'])
    if current['sif_sha256'] != provenance['sif_sha256'] or current['case_manifest_sha256'] != provenance['case_manifest_sha256']:
        raise ValueError('SIF or inputs changed since initial submission')
    if Path(current['sif_path']).resolve() != (Path('/home/helu/containers') / (sif_name + '.sif')).resolve():
        raise ValueError('Submitted SIF name differs from verified path')
    gate_dir = results / 'control_gate'
    gate_dir.mkdir(parents=True, exist_ok=True)
    lock = gate_dir / 'dispatch_started.json'
    # Exclusive creation prevents accidental duplicate submissions on a rerun.
    with lock.open('x') as stream:
        json.dump({'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'provenance': str(provenance_path)}, stream)
    analysis_path = repo / 'plot/python/confluent_wet/cw_pmem_pulse_analysis.py'
    spec = importlib.util.spec_from_file_location('pulse_analysis', analysis_path)
    analysis = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analysis)
    manifest = json.loads(manifest_path.read_text())
    controls_manifest = dict(manifest, cases=[r for r in manifest['cases'] if r['kind'] == 'control'])
    control_path = gate_dir / 'controls_manifest.json'
    control_path.write_text(json.dumps(controls_manifest, indent=2) + '\n')
    summary = analysis.aggregate(results, gate_dir / 'summary', control_path, make_figures=False)
    groups, selected = choose_groups(manifest, summary)
    gate = {'groups': groups, 'selected_pulses': len(selected),
            'held_pulses': sum(r['kind'] == 'pulse' for r in manifest['cases']) - len(selected),
            'selected_cases': [r['case'] for r in selected],
            'interpretation': 'Require all four control seeds to pass within each tau_m/initialization group. Held groups need reviewed longer cases generated locally and synced through pheno.'}
    (gate_dir / 'gate.json').write_text(json.dumps(gate, indent=2) + '\n')
    if not selected:
        print(json.dumps({'selected_pulses': 0, 'held_pulses': gate['held_pulses'], 'gate': str(gate_dir / 'gate.json')}))
        return
    selection = gate_dir / 'selected_inputs.txt'
    selection.write_text(''.join(str(manifest_path.parent / r['case'] / 'run.dat') + '\n' for r in selected))
    env = dict(os.environ)
    env.pop('DEPEND', None)
    env.update(CASE_LIST=str(selection), CONCURRENCY='32', ARRAY_CPUS='16', ARRAY_TIME='04:00:00',
               SKIP_PLOTS='0', PLOT_SCRIPT=str(analysis_path), PLOT_HD_ARGS='',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    output = subprocess.check_output(['bash', str(repo / 'scripts_cluster/submit_array.sh'), sif_name,
                                     str(manifest_path.parent / 'pulses')], cwd=repo, env=env, text=True)
    print(output, flush=True)
    pulse_job = job_id(output)
    gate['pulse_array_job_id'] = pulse_job
    (gate_dir / 'gate.json').write_text(json.dumps(gate, indent=2) + '\n')
    env.update(DEPEND='afterany:' + pulse_job, ANALYSIS_CPUS='2', ANALYSIS_TIME='00:30:00')
    output = subprocess.check_output(['bash', str(repo / 'scripts_cluster/submit_analysis.sh'), str(analysis_path),
                                     str(results), str(results / 'summary'), '--summary', '--manifest',
                                     str(manifest_path)], cwd=repo, env=env, text=True)
    print(output, flush=True)
    gate['summary_job_id'] = job_id(output)
    (gate_dir / 'gate.json').write_text(json.dumps(gate, indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('results', type=Path)
    p.add_argument('sif_name')
    p.add_argument('provenance', type=Path)
    a = p.parse_args()
    launch(a.manifest, a.results, a.sif_name, a.provenance)
