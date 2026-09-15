#!/usr/bin/env python3
"""Chi-only phase explorer: dense coarse stream beside sparse native JSON.

Both movies share one encoded clock. Sparse snapshots are held, not interpolated;
their actual source times remain explicit. No simulations or phase refits.
"""
import argparse
import io
import json
from pathlib import Path

import numpy as np

from cw_mc_tm_scan import TC, parameters
from cw_phase_video_board import sha, save_json

FPS = 25


def render(src, out, selection, original):
    import imageio_ffmpeg
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import colormaps
    from PIL import Image
    import cw_common as cw

    data = json.loads(selection.read_text())
    matches = [r for r in data['runs'] if Path(r['case']).name == src.name]
    assert len(matches) == 1
    row = matches[0]
    oldpath = original/'clips'/row['case']/'video.json'
    old = json.loads(oldpath.read_text())
    p = parameters(src)
    for key in ('case', 'mc', 'tm_over_tc', 'chi0', 'replicate', 'seed', 'tail_mean'):
        assert old[key] == row[key], key
    for key in ('mc', 'chi0', 'seed', 'chi_seed', 'pmem', 'tau_chi', 'nsteps'):
        assert p[key] == row[key], key
    assert p['LX'] == p['LY'] == 256 and p['video_stride'] == 8
    assert p['nvideo'] == 337 and p['ninfo'] == 67400
    assert p['chi_freeze_steps'] == row['preparation_steps']
    assert np.isclose(p['tau_m']/TC, row['tm_over_tc'], rtol=5e-6)
    for name in ('parameters.json', 'video_meta.csv'):
        assert sha(src/name) == old['input_sha256'][name], name
    text = '\n'.join(s for s in (src/'video_meta.csv').read_text().splitlines() if not s.startswith('#'))
    stats = np.atleast_1d(np.genfromtxt(io.StringIO(text), delimiter=',', names=True))
    dense_steps = np.asarray(stats['t'], dtype=int)
    assert np.array_equal(dense_steps, np.arange(p['nsteps']//337+1)*337)
    assert np.isfinite(stats['chi_mean']).all()
    assert np.allclose((dense_steps-row['preparation_steps'])/TC, old['trace_times'], atol=1e-10)
    assert np.array_equal(stats['chi_mean'], old['trace_chi'])
    rawpath = src/'video_chi.u8'
    assert rawpath.stat().st_size == len(dense_steps)*32*32
    raw = np.memmap(rawpath, dtype=np.uint8, mode='r', shape=(len(dense_steps), 32, 32))
    snapshot_steps = np.asarray(old['steps'], dtype=int)
    assert snapshot_steps[-1] == p['nsteps']
    assert list(snapshot_steps) == sorted(int(f.stem[5:]) for f in src.glob('frame*.json') if f.stem[5:].isdigit())
    archive = cw.loadarchive(str(src))
    lut = (colormaps['coolwarm'](np.linspace(0, 1, 256))[:, :3]*255).astype(np.uint8)
    images, means, input_hashes, deviations = [], [], {}, []
    for step in snapshot_steps:
        name = f'frame{step}.json'
        digest = sha(src/name)
        assert digest == old['input_sha256'][name], name
        field = np.asarray(archive.extract_and_read(f'frame{step}')['chi']['value'], dtype=float).reshape(256, 256)
        assert np.isfinite(field).all() and -.000001 <= field.min() <= field.max() <= 1.000001
        # Independent correspondence check: full JSON block means reproduce the stream
        # at coincident timestamps to <= one uint8 level (JSON has finite precision).
        if step % 337 == 0 and step <= dense_steps[-1]:
            block = field.reshape(32, 8, 32, 8).mean(axis=(1, 3))
            q = np.clip(np.floor(block*255+.5), 0, 255).astype(int)
            err = int(np.max(np.abs(q-raw[step//337].astype(int))))
            assert err <= 1, (src, step, err)
            deviations.append(err)
        q = np.clip(np.rint(field*255), 0, 255).astype(np.uint8)
        images.append(lut[np.flipud(q.T)].tobytes())
        means.append(float(field.mean()))
        input_hashes[name] = digest
    # Append the exact final JSON instant if it falls between dense stream samples.
    # The low-resolution panel honestly holds its last available stream sample there.
    steps = np.unique(np.r_[dense_steps, p['nsteps']])
    low_idx = np.searchsorted(dense_steps, steps, side='right')-1
    high_idx = np.searchsorted(snapshot_steps, steps, side='right')-1
    assert np.array_equal(np.unique(high_idx), np.arange(len(snapshot_steps)))
    times = (steps-row['preparation_steps'])/TC
    poster_index = int(np.flatnonzero(times >= times[-1]-500)[0])
    n_encoded = len(steps)+FPS-1
    out.mkdir(parents=True, exist_ok=True)
    movies = {}
    for mode, shape in (('low', [32, 32]), ('high', [256, 256])):
        tmp = out/f'{mode}.tmp.mp4'
        writer = imageio_ffmpeg.write_frames(str(tmp), (256, 256), fps=FPS,
            codec='libx264', pix_fmt_out='yuv420p', macro_block_size=2, quality=None,
            output_params=['-crf', '18', '-preset', 'fast', '-threads', '1', '-movflags', '+faststart'])
        writer.send(None)
        try:
            for i in range(len(steps)):
                if mode == 'low':
                    rgb = lut[np.flipud(raw[low_idx[i]].T)]
                    packed = np.repeat(np.repeat(rgb, 8, axis=0), 8, axis=1).tobytes()
                else:
                    packed = images[high_idx[i]]
                writer.send(packed)
                if i == poster_index:
                    Image.frombytes('RGB', (256, 256), packed).save(out/f'{mode}.png')
            for _ in range(FPS-1):
                writer.send(packed)
        finally:
            writer.close()
        reader = imageio_ffmpeg.read_frames(str(tmp)); probe = next(reader); reader.close()
        count, seconds = imageio_ffmpeg.count_frames_and_secs(str(tmp))
        assert count == n_encoded and abs(seconds-n_encoded/FPS) < .05
        assert probe['size'] == (256, 256) and probe['fps'] == FPS
        assert probe['codec'] == 'h264' and probe['pix_fmt'].startswith('yuv420p')
        path = out/f'{mode}.mp4'; tmp.replace(path)
        movies[mode] = dict(native_shape=shape, width=256, height=256,
            source_frames=len(dense_steps) if mode == 'low' else len(snapshot_steps),
            source_step_interval=337 if mode == 'low' else 67400,
            video_sha256=sha(path), video_bytes=path.stat().st_size,
            poster_sha256=sha(out/f'{mode}.png'))
    input_hashes.update({n: sha(src/n) for n in ('parameters.json', 'video_meta.csv', 'video_chi.u8')})
    meta = {**row, 'fps': FPS, 'frames': len(steps), 'encoded_frames': n_encoded,
            'duration': n_encoded/FPS, 'steps': steps.tolist(), 'times': times.tolist(),
            'video_seconds': (np.arange(len(steps))/FPS).tolist(),
            'sample_steps': dict(low=dense_steps[low_idx].tolist(), high=snapshot_steps[high_idx].tolist()),
            'chi': dict(low=stats['chi_mean'][low_idx].tolist(), high=np.asarray(means)[high_idx].tolist()),
            'trace_times': old['trace_times'], 'trace_chi': old['trace_chi'],
            'snapshot_steps': snapshot_steps.tolist(), 'snapshot_means': means,
            'movies': movies, 'input_sha256': input_hashes, 'source_video_record_sha256': sha(oldpath),
            'render_sha256': sha(Path(__file__)), 'poster_index': poster_index,
            'common_snapshot_checks': len(deviations), 'max_block_quantization_error': max(deviations),
            'colormap': 'coolwarm', 'color_limits': [0, 1],
            'timing_note': 'All dense stream frames retained at 25 fps, shared clock for low/high. '
                           'High holds the most recent real JSON snapshot; no temporal interpolation. '
                           'Exact final JSON time appended; low holds its last actual sample there. '
                           'Final screen held for one second; sample_steps give actual field times.'}
    save_json(out/'chi_video.json', meta)
    print(json.dumps(dict(case=row['case'], dense_frames=len(dense_steps), snapshots=len(snapshot_steps),
                          video_MB=sum(m['video_bytes'] for m in movies.values())/1e6)), flush=True)


def render_all(src, out, selection, original, workers):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    data = json.loads(selection.read_text())
    assert 1 <= workers <= 8
    out.mkdir(parents=True, exist_ok=True)
    completed, failed = [], []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(render, src/r['case'], out/'chi_clips'/r['case'], selection, original): r['case']
                for r in data['runs']}
        for job in as_completed(jobs):
            case = jobs[job]
            try:
                job.result(); completed.append(case)
            except Exception as exc:
                failed.append(dict(case=case, error=repr(exc)))
                print(json.dumps(failed[-1]), flush=True)
            save_json(out/'chi_render_status.json', dict(completed=sorted(completed), failed=failed, total=len(jobs)))
    if failed:
        raise RuntimeError(f'{len(failed)} chi renders failed')


def board(out, selection):
    from matplotlib import colormaps
    data = json.loads(selection.read_text())
    runs = []
    for row in data['runs']:
        root = out/'chi_clips'/row['case']
        r = json.loads((root/'chi_video.json').read_text())
        for k in ('case', 'mc', 'tm_over_tc', 'chi0', 'seed', 'replicate', 'representative', 'tail_mean'):
            assert r[k] == row[k], (row['case'], k)
        assert r['steps'][0] == 0 and r['steps'][-1] == r['nsteps']
        assert len(r['steps']) == len(r['times']) == len(r['video_seconds']) == r['frames']
        for mode, native in (('low', [32, 32]), ('high', [256, 256])):
            m = r['movies'][mode]
            assert m['native_shape'] == native and m['width'] == m['height'] == 256
            assert sha(root/f'{mode}.mp4') == m['video_sha256']
            assert sha(root/f'{mode}.png') == m['poster_sha256']
            assert len(r['sample_steps'][mode]) == len(r['chi'][mode]) == r['frames']
            assert all(s <= t for s, t in zip(r['sample_steps'][mode], r['steps']))
            m['url'] = (root/f'{mode}.mp4').relative_to(out).as_posix()
            m['poster'] = (root/f'{mode}.png').relative_to(out).as_posix()
        runs.append(r)
    for rep in data['representatives']:
        rr = [r for r in runs if r['representative'] == rep['id']]
        assert len(rr) == 4
        assert len({(tuple(r['steps']), r['duration'], r['preparation_steps']) for r in rr}) == 1
    data['runs'] = runs
    data['layout'] = 'low initialization 0, low initialization 1 | high initialization 0, high initialization 1'
    data['chi_colors'] = ['#'+''.join(f'{int(c*255):02x}' for c in rgba[:3])
                          for rgba in colormaps['coolwarm'](np.linspace(0, 1, 256))]
    data['selection_sha256'] = sha(selection)
    template = Path(__file__).with_suffix('.html.in')
    body = template.read_text()
    assert body.count('__DATA__') == 1
    (out/'index.html').write_text(body.replace('__DATA__', json.dumps(data, ensure_ascii=False,
                          allow_nan=False).replace('</', '<\\/')))
    save_json(out/'chi_dashboard_manifest.json', {**data, 'runs': [
        {k: v for k, v in r.items() if k not in ('chi', 'trace_times', 'trace_chi', 'sample_steps', 'times', 'video_seconds')} for r in runs],
        'builder_sha256': sha(Path(__file__)), 'template_sha256': sha(template)})
    print(json.dumps(dict(index=str(out/'index.html'), runs=len(runs), videos=2*len(runs),
                         video_MB=sum(m['video_bytes'] for r in runs for m in r['movies'].values())/1e6)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    for stage in ('render', 'render-all', 'board'):
        p = sub.add_parser(stage)
        for key in (('out', 'selection') if stage == 'board' else ('input', 'out', 'selection', 'original')):
            p.add_argument('--'+key, type=Path, required=True)
        if stage == 'render-all':
            p.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    if args.stage == 'board':
        board(args.out, args.selection)
    elif args.stage == 'render':
        render(args.input, args.out, args.selection, args.original)
    else:
        render_all(args.input, args.out, args.selection, args.original, args.workers)
