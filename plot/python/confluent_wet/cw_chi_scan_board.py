#!/usr/bin/env python3
"""Uniform-chi basin scan: render dense stored streams, then assemble a local video board.

render: inputdir outdir --stage render (submit_analysis_array.sh contract).
board: downloaded results tree, local output directory --stage board.
Uses campaign-wide fixed scales and a shared 136-step video cadence. Numerical
captions use full-resolution video_meta.csv, never quantised movie pixels.
"""
import argparse
import json
from pathlib import Path

import numpy as np

TAU_C = 674.3290333006435
SIGMA_P = 0.03367688582428838
U_RMS = 0.0075639355129229775
EVERY, FPS = 4, 25


def render(src, out):
    import imageio_ffmpeg
    from matplotlib import colormaps
    from PIL import Image, ImageDraw
    import cw_common as cw
    from cw_stream import Stream, _fit_font

    out.mkdir(parents=True, exist_ok=True)
    p = cw.read_params(str(src))
    st = Stream(str(src), p)
    assert int(p['nvideo']) == 34 and int(p['nsteps']) == 202300
    assert len(set(np.diff(st.steps))) == 1 and np.diff(st.steps)[0] == 34
    assert st.steps[-1] >= 202266, 'Incomplete simulation stream'
    panels = [('u', 'magma', 0., 3*U_RMS, '|u|'),
              ('P', 'RdBu_r', -2*SIGMA_P, 2*SIGMA_P, 'P'),
              ('m', 'viridis', 0., 1., 'm'), ('chi', 'bwr', 0., 1., 'chi')]
    # All panels are in one file: their clocks cannot drift apart.
    gap, header, footer = 6, 24, 38
    width, height = st.nx*4+gap*3, st.ny+header+footer
    assert width % 2 == 0 and height % 2 == 0
    canvas = Image.new('RGB', (width, height), '#101318')
    d = ImageDraw.Draw(canvas)
    font = _fit_font('0.0000', st.nx, 14)
    luts, raws = {}, {}
    for j,(name,cmap,low,high,label) in enumerate(panels):
        x=j*(st.nx+gap)
        lut=(colormaps[cmap](np.linspace(0,1,256))[:,:3]*255).astype(np.uint8)
        luts[name], raws[name] = lut, st.raw(name)
        d.text((x+st.nx//2, 12), label, font=font, fill='white', anchor='mm')
        ramp=lut[np.rint(np.linspace(0,255,st.nx)).astype(int)]
        canvas.paste(Image.fromarray(np.tile(ramp[None,:,:],(9,1,1))), (x,header+st.ny+5))
        d.text((x, height-12), f'{low:.3g}', font=font, fill='#cbd5e1', anchor='lm')
        d.text((x+st.nx-2, height-12), f'{high:.3g}', font=font, fill='#cbd5e1', anchor='rm')
    background=np.asarray(canvas).copy()
    selected=np.arange(0,st.n,EVERY)
    path=out/'fields.mp4'
    writer=imageio_ffmpeg.write_frames(str(path),(width,height),fps=FPS,codec='libx264',
            pix_fmt_out='yuv420p',macro_block_size=2,quality=None,
            output_params=['-crf','27','-preset','fast','-threads','2','-movflags','+faststart'])
    writer.send(None)
    try:
        for i in selected:
            frame=background.copy()
            for j,(name,cmap,low,high,label) in enumerate(panels):
                slo,shi=st.limits(name)
                values=slo+raws[name][i].astype(np.float32)*(shi-slo)/255
                q=np.clip(np.rint((values-low)/(high-low)*255),0,255).astype(np.uint8)
                x=j*(st.nx+gap)
                frame[header:header+st.ny,x:x+st.nx]=luts[name][np.flipud(q.T)]
            writer.send(frame.tobytes())
    finally:
        writer.close()
    chi=np.asarray(st.meta['chi_mean'][:st.n])
    window=chi[st.steps >= st.steps[-1]-100*TAU_C]
    half=len(window)//2
    info=dict(case=src.name,chi0=float(p['chi0']),m0=float(p['m0']),
        tau_m=float(p['tau_m'])/TAU_C,tau_c=TAU_C,fps=FPS,frame_steps=136,
        frames=len(selected),duration=len(selected)/FPS,
        times=(st.steps[selected]/TAU_C).round(6).tolist(),
        chi=chi[selected].round(6).tolist(),
        mean_tail=float(window.mean()),std_tail=float(window.std()),
        drift=float(window[half:].mean()-window[:half].mean()),
        last=float(chi[-1]),clip_fraction=st.clip_fraction())
    (out/'video.json').write_text(json.dumps(info))
    print(f'{src.name}: {len(selected)} frames, {path.stat().st_size/1e6:.1f} MB, tail chi={window.mean():.5f}',flush=True)


def board(src,out):
    import shutil
    runs=[]
    out.mkdir(parents=True,exist_ok=True)
    for meta in sorted(src.glob('*/dashboard/video.json')):
        r=json.loads(meta.read_text())
        clip=meta.with_name('fields.mp4')
        if not clip.is_file():
            raise RuntimeError(f'Missing {clip}')
        dest=out/'clips'/f'{r["case"]}.mp4'
        dest.parent.mkdir(exist_ok=True)
        if clip.resolve()!=dest.resolve(): shutil.copy2(clip,dest)
        r['url']=f'clips/{dest.name}'
        runs.append(r)
    assert len(runs)==33, f'Expected 33 complete movies, got {len(runs)}'
    assert len({(r['frames'],r['fps'],r['frame_steps']) for r in runs})==1
    for tm in (12,20,30):
        assert sorted(round(r['chi0'],1) for r in runs if round(r['tau_m'])==tm)==[i/10 for i in range(11)]
    html=Path(__file__).with_suffix('.html').read_text().replace('__DATA__',json.dumps(runs))
    (out/'index.html').write_text(html)
    (out/'summary.json').write_text(json.dumps([{k:v for k,v in r.items() if k not in ('times','chi')} for r in runs],indent=2))
    print(out/'index.html')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('inputdir',type=Path)
    ap.add_argument('outdir',type=Path)
    ap.add_argument('--stage',choices=['render','board'],default='render')
    a=ap.parse_args()
    (render if a.stage=='render' else board)(a.inputdir,a.outdir)
