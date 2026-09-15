#!/usr/bin/env python3
"""Build an offline experimental phase map linked to representative run videos.

Selection and original classifications are frozen from the existing comparison.
Rendering reads every full-resolution JSON snapshot, including the exact final step. No simulation,
phase classification, calibration or theoretical boundary is recomputed.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from cw_mc_tm_scan import TC, parameters

REPRESENTATIVES = [
    ('A', .14, 16., '高 activity', '两种初始化都落入低 χ 状态'),
    ('B', .21, 3., '中间态', '两种初始化都收敛到 χ ≈ 0.5'),
    ('C', .30, 16., '低 activity', '两种初始化都落入高 χ 状态'),
    ('D', .22, 20., 'History dependence', '两种初始化选择不同的长期状态'),
    ('E', .22, 8., '转变附近 · 单态侧', '两种初始化汇合到同一状态'),
    ('F', .22, 12., '转变附近 · 初值相关侧', '相同 mc，较长 memory 下保留两支状态'),
]
PANELS = [('u', 'magma', 0., 3*.0075639355129229775, '|u|'),
          ('P', 'RdBu_r', -2*.03367688582428838, 2*.03367688582428838, 'P'),
          ('m', 'viridis', 0., 1., 'm'), ('chi', 'coolwarm', 0., 1., 'χ')]
FPS = 25


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def select(summary_path, comparison_path, manifest_path, out):
    summary, comparison, manifest = [json.loads(p.read_text())
                                    for p in (summary_path, comparison_path, manifest_path)]
    if summary['missing'] or summary['invalid']:
        raise ValueError('The original phase scan is incomplete')
    groups = [r for r in comparison['predictions']['gaussian_pressure']
              if r['L'] == 256 and 1 <= r['tm_over_tc'] <= 30]
    assert len(groups) == 154
    originals = {r['case']: r for r in manifest['cases']}
    reps, runs = [], []
    for label, mc, tm, title, description in REPRESENTATIVES:
        match = [g for g in groups if g['mc'] == mc and g['tm_over_tc'] == tm]
        assert len(match) == 1
        rows = [r for r in summary['cases']
                if r['L'] == 256 and r['mc'] == mc and r['tm_over_tc'] == tm]
        assert len(rows) == 4
        assert {(r['chi0'], r['replicate']) for r in rows} == {(c, s) for c in (0, 1) for s in (1, 2)}
        assert all(r['complete'] and r['settled_by_diagnostics'] and not r['diagnostic_flags'] for r in rows)
        reps.append(dict(id=label, mc=mc, tm_over_tc=tm, title=title, description=description,
                         classification=match[0]['classification'], unique_chi=match[0]['unique_chi']))
        for row in rows:
            source = originals[row['case']]
            card = manifest_path.parent / row['case'] / 'run.dat'
            assert sha(card) == source['run_dat_sha256']
            pars = dict(tuple(v.strip() for v in line.split('=', 1))
                        for line in card.read_text().splitlines()
                        if line.strip() and not line.lstrip().startswith('#'))
            assert pars['video-stride'] == '8' and pars['nvideo'] == '337'
            runs.append({**row, 'representative': label,
                         'chi_seed': int(pars['chi-seed']),
                         'source_runcard_sha256': source['run_dat_sha256']})
    assert len(runs) == 24
    data = dict(campaign='20260909/cw_mc_tm_phase', tau_c=TC, representatives=reps,
                groups=groups, runs=runs, curve={k: comparison['curves']['gaussian_pressure'][k] for k in ('g', 'folds')},
                palette=dict(active='#63B8C6', middle='#F5F3F0', passive='#E8A0A6',
                             history='#8C6BB1', boundary='#30363B'),
                sources={str(p.resolve()): sha(p) for p in (summary_path, comparison_path, manifest_path)},
                interpretation='Original finite-time phase map; seed_sensitive is purple for display only. '
                               'High/intermediate/low describe the continuous chi scale, not new phase labels.')
    out.mkdir(parents=True, exist_ok=True)
    save_json(out/'selection.json', data)
    files = [f'{r["case"]}/{name}' for r in runs
             for name in ['parameters.json', 'video_meta.csv'] +
             [f'frame{step}.json' for step in sorted(set(range(0, r['nsteps']+1, 67400)) | {r['nsteps']})]]
    (out/'input_filelist.txt').write_text('\n'.join(files)+'\n')
    print(json.dumps({'representatives': len(reps), 'cases': len(runs), 'selection': str(out/'selection.json')}))


def render(src, out, selection):
    import io
    import imageio_ffmpeg
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import colormaps
    from PIL import Image, ImageDraw, ImageFont
    from matplotlib.font_manager import findfont, FontProperties
    import cw_common as cw
    from cw_dash import seed_grid

    data = json.loads(selection.read_text())
    rows = [r for r in data['runs'] if Path(r['case']).name == src.name]
    if len(rows) != 1:
        raise ValueError('Input is not a unique selected original run')
    row, p = rows[0], parameters(src)
    for key in ('mc', 'seed', 'chi_seed', 'chi0', 'nsteps', 'pmem', 'tau_chi'):
        if p[key] != row[key]:
            raise ValueError(f'Unexpected {key}: {src}')
    assert p['LX'] == p['LY'] == 256 and p['chi_config'] == 'uniform'
    assert p['chi_freeze_steps'] == row['preparation_steps']
    assert np.isclose(p['tau_m']/TC, row['tm_over_tc'], rtol=5e-6)
    assert p['ninfo'] == 67400 and p['nstart'] == 0
    steps = sorted(int(f.stem[5:]) for f in src.glob('frame*.json') if f.stem[5:].isdigit())
    expected = sorted(set(range(0, int(p['nsteps'])+1, 67400)) | {int(p['nsteps'])})
    assert steps == expected, 'Missing or extra full-resolution snapshots'
    times = (np.asarray(steps)-row['preparation_steps'])/TC
    text = '\n'.join(s for s in (src/'video_meta.csv').read_text().splitlines() if not s.startswith('#'))
    stats = np.atleast_1d(np.genfromtxt(io.StringIO(text), delimiter=',', names=True))
    assert np.array_equal(stats['t'], np.arange(int(p['nsteps'])//337+1)*337)
    assert all(np.isfinite(stats[name]).all() for name in stats.dtype.names)
    trace_times = (stats['t']-row['preparation_steps'])/TC
    end = (row['nsteps']-row['preparation_steps'])/TC
    tail = trace_times >= end-max(500., 20*row['tm_over_tc'])
    assert abs(float(np.mean(stats['chi_mean'][tail]))-row['tail_mean']) < 1e-9

    # One stored ~100-tau_c interval per second; hold pictures, never interpolate fields.
    # Quantizing each interval to >=1 encoded frame retains close final snapshots too.
    holds = np.r_[np.maximum(1, np.rint(np.diff(steps)/TC/100*FPS).astype(int)), FPS]
    video_seconds = np.r_[0, np.cumsum(holds[:-1])]/FPS
    n_encoded = int(np.sum(holds))
    out.mkdir(parents=True, exist_ok=True)
    panel, gap, header, footer = 512, 16, 40, 52
    width, height = 4*panel+3*gap, panel+header+footer
    canvas = Image.new('RGB', (width, height), '#111820')
    draw, luts = ImageDraw.Draw(canvas), {}
    font = ImageFont.truetype(findfont(FontProperties(family='DejaVu Sans')), 22)
    for j, (name, cmap, low, high, label) in enumerate(PANELS):
        x = j*(panel+gap)
        lut = (colormaps[cmap](np.linspace(0, 1, 256))[:, :3]*255).astype(np.uint8)
        luts[name] = lut
        draw.text((x+panel/2, header/2), label, font=font, anchor='mm', fill='white')
        ramp = lut[np.rint(np.linspace(0, 255, panel)).astype(int)]
        canvas.paste(Image.fromarray(np.tile(ramp[None, :, :], (12, 1, 1))), (x, header+panel+7))
        draw.text((x+2, height-14), f'{low:.3g}', font=font, anchor='lm', fill='#cad4df')
        draw.text((x+panel-2, height-14), f'{high:.3g}', font=font, anchor='rm', fill='#cad4df')
    background = np.asarray(canvas).copy()
    fig = plt.figure(figsize=(panel/100, panel/100), dpi=100)
    fig.patch.set_alpha(0)
    ax = fig.add_axes([0, 0, 1, 1]); ax.patch.set_alpha(0)
    ax.set(xlim=(-.5, 255.5), ylim=(-.5, 255.5)); ax.set_axis_off()
    seeds = seed_grid(256, 24.)
    archive = cw.loadarchive(str(src))
    dest, tmp = out/'fields.mp4', out/'fields.tmp.mp4'
    writer = imageio_ffmpeg.write_frames(str(tmp), (width, height), fps=FPS,
                codec='libx264', pix_fmt_out='yuv420p', macro_block_size=2, quality=None,
                output_params=['-crf', '18', '-preset', 'fast', '-threads', '2', '-movflags', '+faststart'])
    writer.send(None)
    chi, memory, input_hashes = [], [], {}
    poster_index = int(np.flatnonzero(times >= end-500)[0])
    try:
        for i, step in enumerate(steps):
            # read_frame(i) would miss the final snapshot when nsteps % ninfo != 0.
            raw = archive.extract_and_read('frame'+str(step))
            fields = {}
            for name, key in (('ux', 'ux_mat'), ('uy', 'uy_mat'), ('P', 'pressure'), ('m', 'm'), ('chi', 'chi')):
                values = np.asarray(raw[key]['value'], dtype=float)
                assert values.size == 256*256 and np.isfinite(values).all(), (step, name)
                fields[name] = values.reshape(256, 256)
            fields['u'] = np.hypot(fields['ux'], fields['uy'])
            chi.append(float(fields['chi'].mean())); memory.append(float(fields['m'].mean()))
            assert -.000001 <= fields['chi'].min() <= fields['chi'].max() <= 1.000001
            frame = background.copy()
            for j, (name, _, low, high, _) in enumerate(PANELS):
                q = np.clip(np.rint((fields[name]-low)/(high-low)*255), 0, 255).astype(np.uint8)
                rgb = np.repeat(np.repeat(luts[name][np.flipud(q.T)], 2, axis=0), 2, axis=1)
                if name == 'u':
                    for artist in list(ax.collections)+list(ax.patches):
                        artist.remove()
                    assert not ax.collections and not ax.patches
                    if fields['u'].max() > 0:
                        ax.streamplot(np.arange(256), np.arange(256), fields['ux'].T, fields['uy'].T,
                                      start_points=seeds, density=1.4, linewidth=.7, arrowsize=.6, color='white')
                    fig.canvas.draw()
                    overlay = Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy())
                    rgb = np.asarray(Image.alpha_composite(Image.fromarray(rgb).convert('RGBA'), overlay).convert('RGB'))
                x = j*(panel+gap)
                frame[header:header+panel, x:x+panel] = rgb
            if i == poster_index:
                Image.fromarray(frame).save(out/'poster.png')
            packed = frame.tobytes()
            for _ in range(int(holds[i])):
                writer.send(packed)
            input_hashes[f'frame{step}.json'] = sha(src/f'frame{step}.json')
    finally:
        writer.close(); plt.close(fig)
    reader = imageio_ffmpeg.read_frames(str(tmp)); probe = next(reader); reader.close()
    n, seconds = imageio_ffmpeg.count_frames_and_secs(str(tmp))
    assert n == n_encoded and probe['size'] == (width, height) and probe['fps'] == FPS
    assert probe['codec'] == 'h264' and probe['pix_fmt'].startswith('yuv420p')
    assert abs(seconds-n_encoded/FPS) < .05
    tmp.replace(dest)
    input_hashes.update({name: sha(src/name) for name in ('parameters.json', 'video_meta.csv')})
    meta = {**row, 'fps': FPS, 'frames': len(steps), 'encoded_frames': n_encoded,
            'duration': n_encoded/FPS, 'frame_steps': 67400, 'steps': steps,
            'times': times.tolist(), 'video_seconds': video_seconds.tolist(),
            'native_shape': [256, 256], 'block_size': 1, 'source_kind': 'full_resolution_json',
            'width': width, 'height': height, 'chi': chi, 'memory': memory,
            'trace_times': trace_times.tolist(), 'trace_chi': stats['chi_mean'].tolist(),
            'video_sha256': sha(dest), 'video_bytes': dest.stat().st_size,
            'poster_sha256': sha(out/'poster.png'), 'poster_index': poster_index,
            'input_sha256': input_hashes, 'render_sha256': sha(Path(__file__)),
            'panels': PANELS, 'streamlines': {'seed_step': 24, 'density': 1.4},
            'timing_note': 'Each real JSON snapshot is held until the next; nominal 100 tau_c/second. '
                           'Durations are rounded to 1/25 second with at least one encoded frame per snapshot. '
                           'The exact final snapshot is held for one second. No field interpolation.'}
    save_json(out/'video.json', meta)
    print(json.dumps({'case': row['case'], 'snapshots': len(steps), 'encoded_frames': n_encoded,
                      'video_MB': meta['video_bytes']/1e6}), flush=True)


def board(out, selection):
    data = json.loads(selection.read_text())
    runs = []
    for row in data['runs']:
        root = out/'clips'/row['case']
        r = json.loads((root/'video.json').read_text())
        for key in ('case', 'mc', 'tm_over_tc', 'chi0', 'seed', 'replicate', 'representative', 'tail_mean'):
            assert r[key] == row[key], (row['case'], key)
        assert sha(root/'fields.mp4') == r['video_sha256']
        assert (root/'fields.mp4').stat().st_size == r['video_bytes']
        assert r['frames'] == len(r['steps']) == len(r['times']) == len(r['video_seconds']) == len(r['chi']) == len(r['memory'])
        assert r['steps'][0] == 0 and r['steps'][-1] == r['nsteps']
        assert r['native_shape'] == [256, 256] and r['source_kind'] == 'full_resolution_json'
        assert sha(root/'poster.png') == r['poster_sha256']
        r['url'] = (root/'fields.mp4').relative_to(out).as_posix()
        r['poster'] = (root/'poster.png').relative_to(out).as_posix()
        runs.append(r)
    for rep in data['representatives']:
        rr = [r for r in runs if r['representative'] == rep['id']]
        assert len(rr) == 4
        assert len({(tuple(r['steps']), tuple(r['video_seconds']), r['preparation_steps'], r['fps']) for r in rr}) == 1
    data['runs'] = runs
    data['selection_sha256'] = sha(selection)
    template = Path(__file__).with_suffix('.html.in').read_text()
    assert template.count('__DATA__') == 1
    payload = json.dumps(data, ensure_ascii=False, allow_nan=False).replace('</', '<\\/')
    (out/'index.html').write_text(template.replace('__DATA__', payload))
    save_json(out/'dashboard_manifest.json', {**data, 'runs': [
        {k: v for k, v in r.items() if k not in ('chi', 'memory', 'trace_times', 'trace_chi')} for r in runs],
        'builder_sha256': sha(Path(__file__)), 'template_sha256': sha(Path(__file__).with_suffix('.html.in'))})
    print(json.dumps({'index': str(out/'index.html'), 'videos': len(runs),
                      'video_MB': sum(r['video_bytes'] for r in runs)/1e6}))


def render_all(src, out, selection, workers):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    data = json.loads(selection.read_text())
    assert 1 <= workers <= 8
    out.mkdir(parents=True, exist_ok=True)
    completed, failed = [], []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(render, src/r['case'], out/'clips'/r['case'], selection): r['case']
                for r in data['runs']}
        for job in as_completed(jobs):
            case = jobs[job]
            try:
                job.result()
                completed.append(case)
            except Exception as exc:
                failed.append({'case': case, 'error': repr(exc)})
                print(json.dumps(failed[-1]), flush=True)
            save_json(out/'render_status.json', dict(completed=sorted(completed), failed=failed,
                      total=len(jobs), source_kind='full_resolution_json'))
    if failed:
        raise RuntimeError(f'{len(failed)} video renders failed; inspect render_status.json')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='stage', required=True)
    s = sub.add_parser('select')
    for name in ('summary', 'comparison', 'manifest', 'out'):
        s.add_argument('--'+name, type=Path, required=True)
    r = sub.add_parser('render')
    for name in ('input', 'out', 'selection'):
        r.add_argument('--'+name, type=Path, required=True)
    b = sub.add_parser('board')
    for name in ('out', 'selection'):
        b.add_argument('--'+name, type=Path, required=True)
    batch = sub.add_parser('render-all')
    for name in ('input', 'out', 'selection'):
        batch.add_argument('--'+name, type=Path, required=True)
    batch.add_argument('--workers', type=int, default=8)
    a = p.parse_args()
    if a.stage == 'select':
        select(a.summary, a.comparison, a.manifest, a.out)
    elif a.stage == 'render':
        render(a.input, a.out, a.selection)
    elif a.stage == 'render-all':
        render_all(a.input, a.out, a.selection, a.workers)
    else:
        board(a.out, a.selection)
