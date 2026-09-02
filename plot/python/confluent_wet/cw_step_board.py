#!/usr/bin/env python3
"""One dashboard for the whole tau_m axis at the symmetric operating point.

    cw_step_board.py <fine_results> <sym_results> --videos DIR --out DIR

Both trees run at mc = 0.2100, so 0.3-10 from the fine scan and 15/22/30 from the symmetry
test are one continuous axis and are presented as one.

Deliberately spare: a row of tau_m buttons, an `analysis` button, and nothing else. The
selected tau_m shows its clips side by side, labelled only by their initial condition. The
analysis view holds <chi> and std(chi) against tau_m.
"""
import argparse
import base64
import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

ARMS = [("a", "chi = 0"), ("b", "chi = 1"), ("c1", "chi = noise")]


def load(d, pat, clip):
    """{(tau_m, arm): (clip_basename, part)} for one results tree."""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        c = os.path.basename(os.path.dirname(p))
        m = re.match(pat, c)
        if not m:
            continue
        g, arm = float(m.group(1).replace("p", ".")), m.group(2)
        with open(p) as fh:
            out[(g, arm)] = (clip.format(case=c), json.load(fh))
    return out


def figures(R, gs, outdir):
    """<chi> and std(chi) against tau_m, as two standalone PNGs."""
    made = {}
    for key, ylab, fname, logy in (
            ("chi", r"$\langle\chi\rangle$", "chi_vs_taum.png", False),
            ("std", r"$\mathrm{std}(\chi)$", "stdchi_vs_taum.png", False)):
        fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
        for arm, lab in ARMS:
            x = [g for g in gs if (g, arm) in R]
            if not x:
                continue
            y = [R[(g, arm)][1] if key == "chi" else R[(g, arm)][2] for g in x]
            ax.plot(x, y, "-o", lw=1.4, ms=4, label=lab)
        ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel(r"$\tau_m/\tau_c$")
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)
        fig.savefig(os.path.join(outdir, fname), dpi=200)
        plt.close(fig)
        with open(os.path.join(outdir, fname), "rb") as fh:
            made[key] = base64.b64encode(fh.read()).decode()
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fine")
    ap.add_argument("sym")
    ap.add_argument("--videos", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    fine = load(a.fine, r"tm(.+)_([abc]\d?)$", "b1_{case}_chi.mp4")
    sym = load(a.sym, r"s1_tm(.+)_([ab])$", "sym_{case}_chi.mp4")

    # fine wins where they overlap (tau_m = 10): it carries the noise start as well
    R = {}
    for (g, arm), (clip, d) in list(sym.items()) + list(fine.items()):
        s = d.get("settled") or {}
        f = d.get("fate") or {}
        R[(g, arm)] = (clip, f.get("chi_tail", s.get("chi_mean_tail", float("nan"))),
                       d["flow"]["std_chi"])
    gs = sorted(set(g for g, _ in R))

    imgs = figures(R, gs, a.out)

    have = set(os.listdir(a.videos)) if os.path.isdir(a.videos) else set()
    panes = []
    for g in gs:
        cells = []
        for arm, lab in ARMS:
            hit = R.get((g, arm))
            if not hit or hit[0] not in have:
                continue
            cells.append(f'<figure><video data-src="clips/{hit[0]}" muted playsinline loop '
                         f'preload="none"></video><figcaption>{lab}</figcaption></figure>')
        if cells:
            panes.append(f'<div class="pane" data-g="{g:g}">'
                         f'<div class="pg">&tau;_m/&tau;_c = {g:g}</div>'
                         f'{"".join(cells)}</div>')

    buttons = "".join(f'<button class="tb" data-g="{g:g}">{g:g}</button>' for g in gs)

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>tau_m board</title>
<style>
:root{{--bg:#0e1014;--fg:#e9ebef;--dim:#7f8796;--line:#262b34;--on:#5c9bd6;--panel:#171a21}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.5 ui-monospace,"SF Mono",Menlo,monospace}}
.wrap{{max-width:1500px;margin:0 auto;padding:18px 20px 60px;display:flex;
  flex-direction:column;gap:14px}}
.bar{{display:flex;flex-wrap:wrap;gap:6px;align-items:center}}
.lab{{color:var(--dim);margin-right:4px}}
button{{background:var(--panel);color:var(--fg);border:1px solid var(--line);border-radius:3px;
  padding:6px 10px;font:inherit;cursor:pointer}}
button:hover{{border-color:var(--on)}}
button[data-on]{{border-color:var(--on);color:var(--on)}}
.sep{{flex:1}}
.pane{{display:none;grid-template-columns:repeat(3,1fr);gap:8px 12px}}
.pane[data-on]{{display:grid}}
.pg{{grid-column:1/-1;color:var(--dim);border-top:1px solid var(--line);padding-top:8px}}
.panes{{display:flex;flex-direction:column;gap:14px}}
@media(max-width:900px){{.pane[data-on]{{grid-template-columns:1fr}}}}
figure{{margin:0;display:flex;flex-direction:column;gap:6px}}
video{{width:100%;display:block;background:#000;aspect-ratio:400/422;object-fit:contain;
  border:1px solid var(--line);border-radius:3px}}
figcaption{{color:var(--dim);text-align:center}}
.ctl{{display:flex;gap:10px;align-items:center;flex-wrap:wrap}}
input[type=range]{{flex:1;min-width:220px;accent-color:var(--on)}}
#an{{display:none;gap:16px;flex-wrap:wrap;justify-content:center}}
#an[data-on]{{display:flex}}
#an img{{max-width:580px;width:100%;height:auto;background:#fff;border-radius:3px}}
</style></head><body><div class="wrap">

<div class="bar"><span class="lab">videos&nbsp; &tau;_m/&tau;_c:</span>{buttons}
  <span class="sep"></span><button id="anb">analysis</button></div>

<div class="ctl"><button id="pp">play</button>
  <input id="sk" type="range" min="0" max="1000" value="0">
  <span id="tt" class="lab">0.0s</span>
  <span class="lab">speed</span>
  <button class="sp" data-r="0.5">0.5&times;</button>
  <button class="sp" data-r="1" data-on>1&times;</button>
  <button class="sp" data-r="2">2&times;</button>
  <button class="sp" data-r="4">4&times;</button></div>

<div class="panes">{"".join(panes)}</div>

<div id="an"><img alt="mean chi against tau_m" src="data:image/png;base64,{imgs['chi']}">
<img alt="std chi against tau_m" src="data:image/png;base64,{imgs['std']}"></div>

</div><script>
const panes=[...document.querySelectorAll('.pane')],tb=[...document.querySelectorAll('.tb')];
const an=document.getElementById('an'),anb=document.getElementById('anb');
const pp=document.getElementById('pp'),sk=document.getElementById('sk'),tt=document.getElementById('tt');
let rate=1,playing=false;

// every open pane contributes; the first one is the clock everything else is held to
function vids(){{return [...document.querySelectorAll('.pane[data-on] video')]}}
function attach(v,t){{
  if(!v.src&&v.dataset.src){{v.src=v.dataset.src;v.preload='auto';v.load();}}
  v.playbackRate=rate;
  // only the seek needs an event, and it needs metadata rather than data. STARTING is left
  // to the watchdog below: `canplay` fires once, so a listener added after it has already
  // gone by never runs -- which is exactly what happens when a pane is opened (fires) and
  // play is pressed a moment later (too late).
  const go=()=>{{if(v.duration)v.currentTime=t;}};
  if(v.readyState>=1)go(); else v.addEventListener('loadedmetadata',go,{{once:true}});
}}
function toggle(g){{
  an.removeAttribute('data-on'); anb.removeAttribute('data-on');
  const p=panes.find(x=>x.dataset.g===g), b=tb.find(x=>x.dataset.g===g);
  if(!p)return;
  if(p.hasAttribute('data-on')){{
    [...p.querySelectorAll('video')].forEach(v=>v.pause());
    p.removeAttribute('data-on'); b.removeAttribute('data-on');
    if(!vids().length){{playing=false; pp.textContent='play';}}
  }} else {{
    const lead=vids()[0];
    const t=(lead&&lead.duration)?lead.currentTime:0;
    p.setAttribute('data-on',''); b.setAttribute('data-on','');
    [...p.querySelectorAll('video')].forEach(v=>attach(v,t));
  }}
}}
tb.forEach(b=>b.addEventListener('click',()=>toggle(b.dataset.g)));
anb.addEventListener('click',()=>{{
  vids().forEach(v=>v.pause()); playing=false; pp.textContent='play';
  panes.forEach(p=>p.removeAttribute('data-on')); tb.forEach(b=>b.removeAttribute('data-on'));
  an.toggleAttribute('data-on'); anb.toggleAttribute('data-on');
}});
pp.addEventListener('click',()=>{{
  const vs=vids(); if(!vs.length)return;
  playing=!playing; pp.textContent=playing?'pause':'play';
  if(!playing)vs.forEach(v=>v.pause());   // starting is the watchdog's job
}});
sk.addEventListener('input',()=>{{
  const f=sk.value/1000;
  vids().forEach(v=>{{if(v.duration)v.currentTime=f*v.duration}});
}});
document.querySelectorAll('.sp').forEach(b=>b.addEventListener('click',()=>{{
  rate=parseFloat(b.dataset.r);
  document.querySelectorAll('.sp').forEach(x=>x.toggleAttribute('data-on',x===b));
  vids().forEach(v=>v.playbackRate=rate);
}}));
setInterval(()=>{{
  const vs=vids();
  // start anything that should be running and is not: a clip still buffering when play was
  // pressed, or one whose pane was opened afterwards. Self-healing, so no event can be missed.
  if(playing)vs.forEach(v=>{{if(v.paused&&v.readyState>=2)v.play().catch(()=>{{}})}});
  const v=vs[0]; if(!v||!v.duration)return;
  sk.value=Math.round(1000*v.currentTime/v.duration);
  tt.textContent=v.currentTime.toFixed(1)+'s';
  // independently decoded clips drift apart over an 89 s loop; pull them back
  vs.slice(1).forEach(o=>{{if(o.duration&&Math.abs(o.currentTime-v.currentTime)>0.25)
    o.currentTime=v.currentTime}});
}},200);
toggle(tb[0].dataset.g);
</script></body></html>
"""
    out = os.path.join(a.out, "index.html")
    with open(out, "w") as fh:
        fh.write(html)
    print(f"wrote {out} ({len(html)/1024:.0f} KB), {len(gs)} tau_m points, "
          f"{sum(len(p.split('<figure')) - 1 for p in panes)} clips referenced")
    for g in gs:
        arms = [arm for arm, _ in ARMS if (g, arm) in R]
        print(f"  {g:6g}  " + " ".join(f"{a_}={R[(g,a_)][1]:.3f}" for a_ in arms))


if __name__ == "__main__":
    main()
