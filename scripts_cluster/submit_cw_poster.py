#!/usr/bin/env python3
"""Verify and submit one 20260917 poster campaign plus its dependent summary.

Call on jed after GitHub pheno sync. The user-authorized production submission
is nonblocking; exclusive records prevent accidental duplicate arrays.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from verify_campaign_reuse import verify

REPO = Path(__file__).resolve().parents[1]
SCRATCH = Path('/scratch/helu/mass_hd')
CONFIG = {
    'distributions': ('20260917/cw_poster_distributions', 'cw_poster_distribution.py', 16, '03:00:00'),
    'pulse10': ('20260917/cw_poster_pulse10', 'cw_poster_pulse.py', 32, '04:00:00'),
    'pulse50': ('20260917/cw_poster_pulse50', 'cw_poster_pulse.py', 32, '04:00:00'),
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def job_id(output):
    values = re.findall(r'^(?:Submitted batch job )?(\d+)(?:;\S+)?$', output, re.M)
    if len(values) != 1:
        raise ValueError('Expected exactly one Slurm job id: '+output)
    return values[0]


def write_record(path, record):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(record, indent=2)+'\n')
    temporary.replace(path)


def launch(kind, expected_commit, previous, sif_name):
    if REPO != Path('/home/helu/mass_hd'):
        raise ValueError('Production launcher must run in /home/helu/mass_hd on jed')
    campaign, analysis_name, concurrency, walltime = CONFIG[kind]
    manifest_path = REPO/'cases'/campaign/'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    raw = SCRATCH/'cases'/campaign
    results = SCRATCH/'results/cases'/campaign
    proof = SCRATCH/'campaign_provenance'/campaign
    if raw.exists() or results.exists() or (proof/'submission_started.json').exists():
        raise FileExistsError('Campaign already has outputs or a submission record; inspect before retrying')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    if commit != expected_commit:
        raise ValueError('Local/cluster commit mismatch')
    runtime = ['scripts_cluster/submit_array.sh', 'scripts_cluster/submit_case.sh',
               'scripts_cluster/submit_analysis.sh', 'scripts_cluster/verify_campaign_reuse.py',
               'scripts_cluster/submit_cw_poster.py', 'scripts/gen_cw_poster_campaigns.py',
               'plot/python/confluent_wet/'+analysis_name]
    if kind.startswith('pulse'):
        runtime.append('plot/python/confluent_wet/cw_pmem_pulse_analysis.py')
    reused = [row for row in manifest['cases'] if 'reuse_from' in row]
    fresh = [row for row in manifest['cases'] if 'reuse_from' not in row]
    if reused:
        if any(row['kind'] != 'control' for row in reused) or any(row['kind'] != 'pulse' for row in fresh):
            raise ValueError('Reuse must contain controls only, and new runs pulses only')
        runtime += ['scripts/gen_cw_poster_pulse50.py', 'scripts_cluster/prepare_cw_pulse_reuse.py']
    # Check analysis/workflow AND committed case contents before any scheduling.
    input_names = [str(manifest_path.relative_to(REPO))] + [
        str((manifest_path.parent/row['case']/'run.dat').relative_to(REPO)) for row in manifest['cases']]
    if reused:
        input_names += ['cases/'+manifest['source_campaign']+'/manifest.json']
        input_names += ['cases/'+row['reuse_from']+'/run.dat' for row in reused]
    hashes = {}
    for name in runtime+input_names:
        tracked = subprocess.check_output(['git', 'show', f'HEAD:{name}'], cwd=REPO)
        if (REPO/name).read_bytes() != tracked:
            raise ValueError(f'Uncommitted or changed production input: {name}')
        if name in runtime:
            hashes[name] = digest(tracked)
    sif = Path('/home/helu/containers')/(sif_name+'.sif')
    if not sif.is_file():
        raise FileNotFoundError(sif)
    proof.mkdir(parents=True, exist_ok=True)
    verify(manifest_path, previous, sif, expected_commit, proof/'provenance.json')
    mplconfig = proof/'mplconfig'
    mplconfig.mkdir(exist_ok=True)
    record = dict(campaign=campaign, expected_cases=len(manifest['cases']),
                  new_simulations=len(fresh), reused_controls=len(reused),
                  local_commit=expected_commit, cluster_commit=commit,
                  manifest_sha256=digest(manifest_path.read_bytes()),
                  runtime_sha256=hashes, sif_name=sif_name,
                  raw=str(raw), results=str(results), concurrency=concurrency,
                  cpus_per_task=16, walltime=walltime,
                  status='verified_before_submission')
    with (proof/'submission_started.json').open('x') as stream:
        json.dump(record, stream, indent=2)
    record_path = proof/'submission.json'
    write_record(record_path, record)
    env = dict(os.environ)
    for key in ('DEPEND', 'CASE_LIST'):
        env.pop(key, None)
    env.update(CONCURRENCY=str(concurrency), ARRAY_CPUS='16', ARRAY_TIME=walltime,
               SKIP_PLOTS='0', PLOT_SCRIPT=str(REPO/'plot/python/confluent_wet'/analysis_name),
               PLOT_HD_ARGS='--manifest '+str(manifest_path),
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               MPLCONFIGDIR=str(mplconfig), PYTHONDONTWRITEBYTECODE='1')
    if reused:
        env.update(ANALYSIS_CPUS='4', ANALYSIS_TIME='01:00:00')
        command = ['bash', str(REPO/'scripts_cluster/submit_analysis.sh'),
                   str(REPO/'scripts_cluster/prepare_cw_pulse_reuse.py'),
                   str(manifest_path), '--poster']
        output = subprocess.check_output(command, cwd=REPO, env=env, text=True)
        print(output, flush=True)
        record.update(control_reuse_job_id=job_id(output), control_reuse_command=command,
                      status='control_reuse_submitted')
        write_record(record_path, record)
        selected = proof/'pulse_inputs.txt'
        selected.write_text(''.join(str(manifest_path.parent/row['case']/'run.dat')+'\n' for row in fresh))
        env.update(CASE_LIST=str(selected), DEPEND='afterok:'+record['control_reuse_job_id'])
    command = ['bash', str(REPO/'scripts_cluster/submit_array.sh'), sif_name,
               'cases/'+campaign]
    record['array_command'] = command
    output = subprocess.check_output(command, cwd=REPO, env=env, text=True)
    print(output, flush=True)
    record.update(array_job_id=job_id(output), array_submission_output=output,
                  status='array_submitted')
    write_record(record_path, record)
    env.update(DEPEND='afterany:'+record['array_job_id'], ANALYSIS_CPUS='4',
               ANALYSIS_TIME='02:00:00')
    command = ['bash', str(REPO/'scripts_cluster/submit_analysis.sh'),
               str(REPO/'plot/python/confluent_wet'/analysis_name),
               str(results), str(results/'summary'), '--manifest', str(manifest_path), '--summary']
    record['summary_command'] = command
    output = subprocess.check_output(command, cwd=REPO, env=env, text=True)
    print(output, flush=True)
    record.update(summary_job_id=job_id(output), summary_submission_output=output,
                  status='array_and_summary_submitted')
    write_record(record_path, record)
    print(json.dumps({k: record[k] for k in ('campaign', 'expected_cases', 'array_job_id',
                                           'summary_job_id', 'status')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('kind', choices=CONFIG)
    p.add_argument('expected_commit')
    p.add_argument('--previous', type=Path,
                   default=SCRATCH/'campaign_provenance/20260915/cw_pressure_compare/provenance.json')
    p.add_argument('--sif-name', default='mass_cw_pmem_pulse_20260913')
    args = p.parse_args()
    launch(args.kind, args.expected_commit, args.previous, args.sif_name)
