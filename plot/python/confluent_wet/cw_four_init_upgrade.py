#!/usr/bin/env python3
"""Render forty native streams and upgrade their entries in the existing board."""
import argparse
import json
from pathlib import Path

from cw_four_init_board import render, digest

GRID = {0.3, 6.821052632, 7.752631579, 8.684210526, 18.0}


def completed(summary):
    data = json.loads(summary.read_text())
    assert data['available_cases'] == data['expected_cases'] == 40
    assert not data['missing'] and not data['invalid']
    assert {r['tm_over_tc'] for r in data['cases']} == GRID
    assert all(r['complete'] for r in data['cases'])
    assert {(r['tm_over_tc'], r['initialization'], r['replicate']) for r in data['cases']} == {
        (t, s, r) for t in GRID for s in ('0', '1', 'stripe01', 'noise') for r in (1, 2)}
    return data


def upgrade(src, existing, out, summary):
    data = completed(summary)
    original_page = existing/'index.html'
    original = json.loads(original_page.read_text().split('<script type="application/json" id="data">')[1].split('</script>')[0])
    originals = {r['case']: r for r in original['runs']}
    assert len(originals) == 160
    replacements = {}
    for row in data['cases']:
        root = src/row['case']/'dashboard'
        r = json.loads((root/'video.json').read_text())
        assert row['case'] in originals
        for key in ('case', 'tm_over_tc', 'initialization', 'replicate', 'seed', 'chi_seed', 'chi0',
                    'nsteps', 'preparation_steps', 'tail_mean', 'tail_std_time', 'tail_half_drift',
                    'release_chi_mean', 'release_chi_std', 'diagnostic_flags'):
            assert r[key] == row[key], (row['case'], key)
        for key in ('tm_over_tc', 'initialization', 'replicate', 'seed', 'chi_seed', 'nsteps', 'preparation_steps', 'fps', 'frame_steps'):
            assert r[key] == originals[row['case']][key], (row['case'], key)
        assert r['native_shape'] == [256, 256] and r['block_size'] == 1
        assert r['width'] == 1060 and r['height'] == 332
        assert r['frames'] == len(r['chi']) == len(r['memory']) == row['nsteps']//337+1
        assert r['source_campaign'] == data['campaign'] and r['summary_sha256'] == digest(summary)
        assert r['render_sha256'] == digest(Path(__file__).with_name('cw_four_init_board.py'))
        assert digest(root/'fields.mp4') == r['video_sha256']
        assert (root/'fields.mp4').stat().st_size == r['video_bytes']
        assert digest(root/'poster.png') == r['poster_sha256']
        r['url'] = (Path('fullres')/row['case']/'dashboard/fields.mp4').as_posix()
        r['poster'] = (Path('fullres')/row['case']/'dashboard/poster.png').as_posix()
        r['replaces_video_sha256'] = originals[row['case']]['video_sha256']
        replacements[row['case']] = r
    runs = [replacements.get(r['case'], r) for r in original['runs']]
    for r in runs:
        if r['case'] not in replacements:
            assert r == originals[r['case']] and (existing/r['url']).is_file()
    for tm in {r['tm_over_tc'] for r in runs}:
        rr = [r for r in runs if r['tm_over_tc'] == tm]
        assert len(rr) == 8 and len({(r['frames'], r['fps'], r['frame_steps'], r['preparation_steps'], tuple(r['native_shape'])) for r in rr}) == 1
    groups = {g['tm_over_tc']: g for g in data['groups']}
    payload = {**original, 'runs': runs, 'groups': [groups.get(g['tm_over_tc'], g) for g in original['groups']],
               'fullres_grid': sorted(GRID), 'fullres_campaign': data['campaign']}
    template = Path(__file__).with_name('cw_four_init_board.html.in')
    html = template.read_text()
    assert html.count('__DATA__') == 1
    out.mkdir(parents=True, exist_ok=True)
    (out/'index.html').write_text(html.replace('__DATA__', json.dumps(payload, ensure_ascii=False,
                            allow_nan=False).replace('</', '<\\/')))
    manifest = dict(runs=[{k: v for k, v in r.items() if k not in ('chi', 'memory')} for r in runs],
        summary_sha256=digest(summary), original_index_sha256=digest(original_page),
        original_manifest_sha256=digest(existing/'dashboard_manifest.json'),
        builder_sha256=digest(Path(__file__)), template_sha256=digest(template),
        fullres_campaign=data['campaign'], fullres_grid=sorted(GRID), replaced_cases=sorted(replacements))
    (out/'dashboard_manifest.json').write_text(json.dumps(manifest, indent=2))
    (out/'filelist.txt').write_text('\n'.join(['index.html','dashboard_manifest.json']+
        [r['url'] for r in runs]+[r['poster'] for r in replacements.values()])+'\n')
    print(json.dumps(dict(index=str(out/'index.html'), native_videos=40, preserved_videos=120,
        new_video_MB=sum(r['video_bytes'] for r in replacements.values())/1e6)), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('input', 'out', 'summary', 'existing'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    upgrade(a.input, a.existing, a.out, a.summary)
