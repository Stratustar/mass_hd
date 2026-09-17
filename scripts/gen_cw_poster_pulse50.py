#!/usr/bin/env python3
"""Create 50% pulses paired to byte-identical completed 10% campaign controls."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from gen_cw_poster_campaigns import REPO, immutable_write, parse_card, source_provenance

CAMPAIGN = '20260917/cw_poster_pulse50'
SOURCE = '20260917/cw_poster_pulse10'


def main():
    source_root, root = REPO/'cases'/SOURCE, REPO/'cases'/CAMPAIGN
    source_path = source_root/'manifest.json'
    old = json.loads(source_path.read_text())
    manifest = deepcopy(old)
    manifest.update(source_provenance())
    assert manifest['simulator_source_sha256'] == old['simulator_source_sha256']
    manifest.update(campaign=CAMPAIGN, source_campaign=SOURCE,
                    source_manifest_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                    relative_threshold_decrease=.5, new_simulations=150, reused_controls=150)
    pulse_value = .5*manifest['fixed']['pmem']
    manifest['fixed']['pmem_pulse_value'] = pulse_value
    for row, original in zip(manifest['cases'], old['cases']):
        source_card = source_root/row['case']/'run.dat'
        raw = source_card.read_text()
        assert hashlib.sha256(raw.encode()).hexdigest() == original['run_dat_sha256']
        if row['kind'] == 'control':
            row['reuse_from'] = SOURCE+'/'+row['case']
            text = raw
        else:
            expected = row['expected_parameters']
            expected['pmem-pulse-value'] = f'{pulse_value:.16g}'
            comments = [CAMPAIGN, '50% sensing-threshold pulse; all other inputs inherited from '+SOURCE,
                        'pc=.016838 -> .008419 -> .016838 for3tau_c; no direct mechanical force.',
                        'Same preparation,2000tau_c feedback wait,1000tau_c recovery and3paired seeds.',
                        'Paired control is reused byte-for-byte; average all3 signed responses before fitting.']
            text = ''.join('# '+line+'\n' for line in comments)+''.join(
                f'{key:24s} = {value}\n' for key, value in expected.items())
        actual, prior = parse_card(text), parse_card(raw)
        changed = {key for key in set(actual)|set(prior) if actual.get(key) != prior.get(key)}
        assert changed == ({'pmem-pulse-value'} if row['kind'] == 'pulse' else set())
        assert actual == row['expected_parameters']
        row['run_dat_sha256'] = hashlib.sha256(text.encode()).hexdigest()
        immutable_write(root/row['case']/'run.dat', text)
    manifest['protocol'] = old['protocol'].replace('down10%', 'down50%')+' All150 original zero-duration controls are reused unchanged.'
    immutable_write(root/'manifest.json', json.dumps(manifest, indent=2)+'\n')
    immutable_write(root/'README.md', '# Dense50% sensing-threshold pulse scan\n\n'
        '150new pulses:25memory times x2uniform histories x3seeds. Reuse150completed controls '
        'from cw_poster_pulse10; their input files remain byte-identical, including inactive '
        'pulse-value settings with zero pulse duration.\n\n'
        'Only the applied pulse threshold changes: pc=.016838 to.008419 for3tau_c. '
        'Core mc=.2287,L256,r=.3,tau_chi=202.3,all seeds,2000tau_c feedback wait and1000tau_c '
        'recovery match the10% campaign. No change to the simulator or fit criteria.\n\n'
        'All3signed pulse-minus-control responses are averaged before fitting, independently '
        'for each initial history. Retain drift, weak-signal and persistent-residual outcomes. '
        'This is a finite-amplitude response protocol.\n')
    print(json.dumps({'campaign':CAMPAIGN,'new_pulses':150,'reused_controls':150,
                      'pulse_pc':pulse_value,'only_changed_parameter':'pmem-pulse-value'}))


if __name__ == '__main__':
    main()
