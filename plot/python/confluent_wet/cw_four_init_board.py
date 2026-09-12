#!/usr/bin/env python3
"""Render the 160-run four-start campaign and build a portable offline video board.

Use submit_analysis_array.sh for render; assemble the board locally after fetching
dashboard/{fields.mp4,video.json}. All stored frames are retained at fixed cadence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from cw_mc_tm_scan import TC, parameters
from cw_stream import Stream, _fit_font

FPS = 25
PANELS = [('u', 'magma', 0., 3*.0075639355129229775, '|u|'),
          ('P', 'RdBu_r', -2*.03367688582428838, 2*.03367688582428838, 'P'),
          ('m', 'viridis', 0., 1., 'm'), ('chi', 'coolwarm', 0., 1., 'chi')]
STARTS = ('0', '1', 'stripe01', 'noise')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def render(src, out, summary_path):
    import imageio_ffmpeg
    from matplotlib import colormaps
    from PIL import Image, ImageDraw
    summary = json.loads(summary_path.read_text())
    rows = [r for r in summary['cases'] if Path(r['case']).name == src.name]
    if len(rows) != 1:
        raise ValueError('Case not uniquely identified in completed summary')
    row = rows[0]
    p = parameters(src)
    for key in ('mc', 'seed', 'chi_seed', 'chi0', 'chi_config', 'nsteps'):
        if row[key] != p[key]:
            raise ValueError(f'Summary/parameters mismatch: {key}')
    if not row['complete'] or p['chi_freeze_steps'] != row['preparation_steps']:
        raise ValueError('Incomplete simulation or wrong preparation')
    if not np.isclose(p['tau_m']/TC, row['tm_over_tc'], rtol=5e-6):
        raise ValueError('Wrong memory coordinate')
    st = Stream(str(src), p)
    expected = int(p['nsteps'])//337+1
    if (st.n != expected or len(st.meta['t']) != expected or st.steps[0] != 0
            or st.nvideo != 337 or st.stride != 8 or (st.nx, st.ny) != (32, 32)
            or not np.all(np.diff(st.steps) == 337)):
        raise ValueError('Unexpected/incomplete video grid or sampling cadence')
    for name, *_ in PANELS:
        if (src/f'video_{name}.u8').stat().st_size != expected*32*32:
            raise ValueError(f'Wrong byte count for {name}')
    out.mkdir(parents=True, exist_ok=True)
    panel, gap, header, footer = 128, 8, 24, 38
    width, height = 4*panel+3*gap, panel+header+footer
    canvas = Image.new('RGB', (width, height), '#101318')
    draw = ImageDraw.Draw(canvas)
    font = _fit_font('0.0000', panel, 13)
    luts, raw = {}, {}
    for j, (name, cmap, low, high, label) in enumerate(PANELS):
        x = j*(panel+gap)
        lut = (colormaps[cmap](np.linspace(0, 1, 256))[:, :3]*255).astype(np.uint8)
        luts[name], raw[name] = lut, st.raw(name)
        draw.text((x+panel//2, 12), label, font=font, fill='white', anchor='mm')
        ramp = lut[np.rint(np.linspace(0, 255, panel)).astype(int)]
        canvas.paste(Image.fromarray(np.tile(ramp[None, :, :], (9, 1, 1))), (x, header+panel+5))
        draw.text((x, height-12), f'{low:.3g}', font=font, fill='#cbd5e1', anchor='lm')
        draw.text((x+panel-2, height-12), f'{high:.3g}', font=font, fill='#cbd5e1', anchor='rm')
    background = np.asarray(canvas).copy()
    path, temp = out/'fields.mp4', out/'fields.tmp.mp4'
    writer = imageio_ffmpeg.write_frames(str(temp), (width, height), fps=FPS,
        codec='libx264', pix_fmt_out='yuv420p', macro_block_size=2, quality=None,
        output_params=['-crf', '23', '-preset', 'fast', '-threads', '2', '-movflags', '+faststart'])
    writer.send(None)
    try:
        for i in range(st.n):
            frame = background.copy()
            for j, (name, _, low, high, _) in enumerate(PANELS):
                slo, shi = st.limits(name)
                values = slo+raw[name][i].astype(np.float32)*(shi-slo)/255
                q = np.clip(np.rint((values-low)/(high-low)*255), 0, 255).astype(np.uint8)
                rgb = luts[name][np.flipud(q.T)]
                # Nearest-neighbour enlargement preserves the true 32x32 block grid.
                rgb = np.repeat(np.repeat(rgb, 4, axis=0), 4, axis=1)
                x = j*(panel+gap)
                frame[header:header+panel, x:x+panel] = rgb
            writer.send(frame.tobytes())
    finally:
        writer.close()
    reader = imageio_ffmpeg.read_frames(str(temp))
    probe = next(reader)
    reader.close()
    decoded_frames, decoded_seconds = imageio_ffmpeg.count_frames_and_secs(str(temp))
    if (decoded_frames != st.n or probe['codec'] != 'h264'
            or not probe['pix_fmt'].startswith('yuv420p') or probe['size'] != (width, height)
            or probe['fps'] != FPS or abs(decoded_seconds-st.n/FPS) > .05):
        raise ValueError(f'Encoded movie mismatch: {probe}')
    temp.replace(path)
    info = {k: row[k] for k in ('case', 'tm_over_tc', 'initialization', 'replicate', 'seed',
        'chi_seed', 'chi0', 'nsteps', 'preparation_steps', 'tail_mean', 'tail_std_time',
        'tail_half_drift', 'diagnostic_flags', 'release_chi_mean', 'release_chi_std')}
    info.update(fps=FPS, frames=st.n, duration=st.n/FPS, frame_steps=337, tau_c=TC,
        start_tc=-row['preparation_steps']/TC,
        end_tc=(float(st.steps[-1])-row['preparation_steps'])/TC,
        width=width, height=height, native_shape=[32, 32], block_size=8,
        chi=np.round(st.meta['chi_mean'], 7).tolist(),
        memory=np.round(st.meta['m_mean'], 7).tolist(),
        clip_fraction=st.clip_fraction(), panels=PANELS,
        video_bytes=path.stat().st_size, video_sha256=digest(path),
        render_sha256=digest(Path(__file__)), analysis_sha256=row['script_sha256'])
    (out/'video.json').write_text(json.dumps(info, allow_nan=False))
    print(json.dumps({'case': row['case'], 'frames': st.n, 'MB': path.stat().st_size/1e6}), flush=True)


def board(src, out, summary_path):
    summary = json.loads(summary_path.read_text())
    if summary['available_cases'] != 160 or summary['missing'] or summary['invalid']:
        raise ValueError('Expected complete 160-case analysis')
    runs = []
    for row in summary['cases']:
        meta = src/row['case']/'dashboard/video.json'
        r = json.loads(meta.read_text())
        for key in ('case', 'initialization', 'replicate', 'tm_over_tc', 'seed', 'chi0', 'preparation_steps'):
            if r[key] != row[key]:
                raise ValueError(f'Wrong {key}: {meta}')
        movie = meta.with_name('fields.mp4')
        if digest(movie) != r['video_sha256'] or movie.stat().st_size != r['video_bytes']:
            raise ValueError(f'Movie does not match metadata: {movie}')
        if r['frames'] != row['nsteps']//337+1 or len(r['chi']) != r['frames'] or len(r['memory']) != r['frames']:
            raise ValueError('Frame or scalar count mismatch')
        r['url'] = str(movie.relative_to(out)) if movie.is_relative_to(out) else ''
        if not r['url']:
            target = out/'clips'/f'{Path(row["case"]).name}.mp4'
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(movie, target)
            r['url'] = target.relative_to(out).as_posix()
        runs.append(r)
    keys = {(r['tm_over_tc'], r['initialization'], r['replicate']) for r in runs}
    if len(keys) != 160 or len({(r['fps'], r['frame_steps'], r['native_shape'][0]) for r in runs}) != 1:
        raise ValueError('Duplicate case or inconsistent movie cadence')
    for tm in {r['tm_over_tc'] for r in runs}:
        rr = [r for r in runs if r['tm_over_tc'] == tm]
        if len(rr) != 8 or len({(r['preparation_steps'], r['frames']) for r in rr}) != 1:
            raise ValueError('A memory-time group cannot be synchronized exactly')
    out.mkdir(parents=True, exist_ok=True)
    template = Path(__file__).with_suffix('.html.in').read_text()
    if template.count('__DATA__') != 1:
        raise ValueError('Invalid template marker')
    payload = dict(runs=runs, groups=summary['groups'], tau_c=TC)
    text = template.replace('__DATA__', json.dumps(payload, ensure_ascii=False).replace('</', '<\\/'))
    if '__DATA__' in text:
        raise ValueError('Unexpanded template')
    (out/'index.html').write_text(text)
    (out/'dashboard_manifest.json').write_text(json.dumps({
        'runs': [{k: v for k, v in r.items() if k not in ('chi', 'memory')} for r in runs],
        'summary_sha256': digest(summary_path), 'builder_sha256': digest(Path(__file__)),
        'template_sha256': digest(Path(__file__).with_suffix('.html.in'))}, indent=2))
    (out/'filelist.txt').write_text('\n'.join(['index.html', 'dashboard_manifest.json']+
        [r['url'] for r in runs])+ '\n')
    print(json.dumps({'board': str(out/'index.html'), 'videos': len(runs),
        'video_MB': sum(r['video_bytes'] for r in runs)/1e6, 'html_MB': len(text.encode())/1e6}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputdir', type=Path)
    p.add_argument('outdir', type=Path)
    p.add_argument('--stage', choices=('render', 'board'), default='render')
    p.add_argument('--summary', required=True, type=Path)
    a = p.parse_args()
    (render if a.stage == 'render' else board)(a.inputdir, a.outdir, a.summary)
