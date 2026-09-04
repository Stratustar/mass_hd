#!/usr/bin/env python3
"""Spare video dashboards for the step-switch campaign. Three modes, one skeleton.

    cw_step_board.py taum <fine> <sym>            --videos DIR --out DIR
    cw_step_board.py tchi <tchi> <fine> <sym>     --videos DIR --out DIR
    cw_step_board.py long <long> [<long2> ...]    --videos DIR --out DIR

All three present the same way: a row of tau_m buttons that toggle (several can be open at
once), an `analysis` button, and one shared transport -- play/pause, scrubber, 0.5-4x --
driving every open clip from a common clock. Nothing else; the labels are the explanation.

What differs is what a pane holds:

  taum   one row, the three initial conditions. mc = 0.2100 across 0.3-30 tau_c.
  tchi   three rows, one per tau_chi rule, each with the three starts. The comparison the
         group exists for reads down a column: same tau_m, same start, different phenotype
         clock.
  long   two rows, chi and m, for the 500000-step runs -- m is on the same page because the
         question there is whether the memory is still moving, not just the phenotype.

Analysis views are the PNGs the per-group analyses already wrote, inlined as data URIs.
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
TCHI_RULES = [(0.3, "fixed", "tau_chi = 0.3 tau_c"),
              (1.0, "fixed", "tau_chi = 1.0 tau_c"),
              (0.5, "prop", "tau_chi = 0.5 tau_m")]


def scan(d, pat, clip):
    """[(tau_m, arm, clip_basename, part)] for one results tree."""
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        c = os.path.basename(os.path.dirname(p))
        m = re.match(pat, c)
        if not m:
            continue
        with open(p) as fh:
            out.append((float(m.group(1).replace("p", ".")), m.group(2),
                        clip.format(case=c), json.load(fh)))
    return out


def b64(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def taum_figs(R, gs, outdir):
    """<chi> and std(chi) against tau_m -- the taum mode writes its own analysis."""
    made = []
    for key, ylab, fname in (("chi", r"$\langle\chi\rangle$", "chi_vs_taum.png"),
                             ("std", r"$\mathrm{std}(\chi)$", "stdchi_vs_taum.png")):
        fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
        for arm, lab in ARMS:
            x = [g for g in gs if (g, arm) in R]
            if not x:
                continue
            ax.plot(x, [R[(g, arm)][key] for g in x], "-o", lw=1.4, ms=4, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel(r"$\tau_m/\tau_c$")
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)
        fig.savefig(os.path.join(outdir, fname), dpi=200)
        plt.close(fig)
        made.append(b64(os.path.join(outdir, fname)))
    return made


def cell(clip, lab, have):
    if clip not in have:
        return ('<figure><div class="ph">-</div>'
                f'<figcaption>{lab}</figcaption></figure>')
    return (f'<figure><video data-src="clips/{clip}" muted playsinline loop '
            f'preload="none"></video><figcaption>{lab}</figcaption></figure>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["taum", "tchi", "long", "fields"])
    ap.add_argument("trees", nargs="*",
                    help="results trees; fields mode needs none (it reads the clip names)")
    ap.add_argument("--videos", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--analysis", nargs="*", default=None,
                    help="PNGs for the analysis view; taum draws its own if omitted")
    ap.add_argument("--title", default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    have = set(os.listdir(a.videos)) if os.path.isdir(a.videos) else set()

    panes, gs, imgs = [], [], []

    if a.mode == "taum":
        fine = scan(a.trees[0], r"tm(.+)_([abc]\d?)$", "b1_{case}_chi.mp4")
        sym = scan(a.trees[1], r"s1_tm(.+)_([ab])$", "sym_{case}_chi.mp4")
        R = {}
        for g, arm, clip, d in sym + fine:          # fine wins the tau_m = 10 overlap
            s, f = d.get("settled") or {}, d.get("fate") or {}
            R[(g, arm)] = {"clip": clip,
                           "chi": f.get("chi_tail", s.get("chi_mean_tail", float("nan"))),
                           "std": d["flow"]["std_chi"]}
        gs = sorted(set(g for g, _ in R))
        imgs = taum_figs(R, gs, a.out)
        for g in gs:
            cells = [cell(R[(g, k)]["clip"], lab, have) for k, lab in ARMS if (g, k) in R]
            if cells:
                panes.append((g, [("", cells)]))

    elif a.mode == "tchi":
        tchi = scan(a.trees[0], r"(?:t1|th|t0p3)_tm(.+)_([abc]\d?)$", "tchi_{case}_chi.mp4")
        fine = scan(a.trees[1], r"tm(.+)_([abc]\d?)$", "b1_{case}_chi.mp4")
        sym = scan(a.trees[2], r"s1_tm(.+)_([ab])$", "sym_{case}_chi.mp4")
        tau_c = None
        R = {}
        for g, arm, clip, d in tchi + fine + sym:
            u = d["params"]
            tau_c = tau_c or (u["tau_m"] / g if g else None)
            # classify by the run's OWN clocks, so the three trees pool without a table
            frac, gchi = u["tau_chi"] / u["tau_m"], u["tau_chi"] / tau_c
            rule = ("prop", 0.5) if abs(frac - 0.5) < 0.02 else \
                   ("fixed", round(gchi, 1)) if min(abs(gchi - 1.0), abs(gchi - 0.3)) < 0.05 \
                   else None
            if rule:
                R[(rule, round(g, 3), arm)] = clip
        # the five tau_m the group was designed on
        gs = sorted(set(g for _, g, _ in R if any(
            abs(g - t) / t < 0.06 for t in (0.3, 3, 10, 15, 30))))
        for g in gs:
            rows = []
            for coef, kind, lab in TCHI_RULES:
                key = ("prop", 0.5) if kind == "prop" else ("fixed", coef)
                cells = [cell(R[(key, g, k)], al, have) for k, al in ARMS
                         if (key, g, k) in R]
                if cells:
                    rows.append((lab, cells))
            if rows:
                panes.append((g, rows))

    elif a.mode == "fields":
        # clips are fld_<tag>_<field>.mp4 and there is no part.json to read: the tags carry
        # the tau_m, so the grid is built from the filenames themselves.
        FLD = [("chi", "chi"), ("m", "m"), ("P", "P"),
               ("p_LB", "p_LB"), ("sigma_B", "sigma_B"), ("Qabs", "|Q|")]
        def _num(t):
            # tags are normally tmX; anything else (an initial-condition name, say) sorts
            # as text and is shown verbatim on its button
            try:
                return (0, float(t[2:].replace("p", ".")) if t.startswith("tm")
                        else float(t.replace("p", ".")))
            except ValueError:
                return (1, 0.0)
        tags = sorted({f.split("_")[1] for f in have if f.startswith("fld_")},
                      key=lambda t: (_num(t), t))
        for t in tags:
            g = t[2:].replace("p", ".") if t.startswith("tm") else t
            rows, cells = [], []
            # only the fields that actually exist: this mode serves both the six-field
            # comparison and smaller ad-hoc sets, and a row of "-" placeholders for fields
            # nobody asked to render is just noise
            for fk, flab in [(k, l) for k, l in FLD
                             if f"fld_{t}_{k}.mp4" in have]:
                cells.append(cell(f"fld_{t}_{fk}.mp4", flab, have))
                if len(cells) == 3:
                    rows.append(("", cells)); cells = []
            if cells:
                rows.append(("", cells))
            if rows:
                panes.append((g, rows))
        gs = [g for g, _ in panes]

    else:  # long
        R = {}
        for t in a.trees:
            for g, arm, clip, d in scan(t, r"tm(.+)_([abc]\d?)$", "long_{case}_chi.mp4"):
                R[(g, arm)] = clip
        gs = sorted(set(g for g, _ in R))
        for g in gs:
            rows = []
            for fld, flab in (("chi", "chi"), ("m", "m")):
                cells = [cell(R[(g, k)].replace("_chi.mp4", f"_{fld}.mp4"), al, have)
                         for k, al in ARMS if (g, k) in R]
                if cells:
                    rows.append((flab, cells))
            if rows:
                panes.append((g, rows))

    if a.analysis:
        imgs = [b64(p) for p in a.analysis if os.path.exists(p)]
    if not panes:
        raise RuntimeError("no panes built -- check the trees and the clips directory")

    pane_html = []
    for g, rows in panes:
        lbl = f"&tau;_m/&tau;_c = {g:g}" if isinstance(g, float) else str(g)
        body = [f'<div class="pg">{lbl}</div>']
        for rlab, cells in rows:
            body.append(f'<div class="prow">'
                        + (f'<div class="rl">{rlab}</div>' if rlab else "")
                        + "".join(cells) + "</div>")
        gid = f"{g:g}" if isinstance(g, float) else str(g)
        pane_html.append(f'<div class="pane" data-g="{gid}">{"".join(body)}</div>')
    buttons = "".join(
        f'<button class="tb" data-g="{g:g}">{g:g}</button>' if isinstance(g, float)
        else f'<button class="tb" data-g="{g}">{g}</button>' for g, _ in panes)
    an_html = "".join(f'<img alt="analysis" src="data:image/png;base64,{i}">' for i in imgs)
    title = a.title or {"taum": "tau_m board", "tchi": "tau_chi board",
                        "long": "500k board", "fields": "fields board"}[a.mode]

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
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
.panes{{display:flex;flex-direction:column;gap:14px}}
.pane{{display:none;flex-direction:column;gap:10px}}
.pane[data-on]{{display:flex}}
.pg{{color:var(--dim);border-top:1px solid var(--line);padding-top:8px}}
.prow{{display:grid;grid-template-columns:repeat(3,1fr);gap:6px 12px}}
.rl{{grid-column:1/-1;color:var(--on);font-size:12.5px}}
@media(max-width:900px){{.prow{{grid-template-columns:1fr}}}}
figure{{margin:0;display:flex;flex-direction:column;gap:6px}}
video,.ph{{width:100%;display:block;background:#000;aspect-ratio:400/422;object-fit:contain;
  border:1px solid var(--line);border-radius:3px}}
.ph{{display:grid;place-items:center;background:var(--panel);color:var(--dim)}}
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

<div class="panes">{"".join(pane_html)}</div>

<div id="an">{an_html}</div>

</div><script>
const panes=[...document.querySelectorAll('.pane')],tb=[...document.querySelectorAll('.tb')];
const an=document.getElementById('an'),anb=document.getElementById('anb');
const pp=document.getElementById('pp'),sk=document.getElementById('sk'),tt=document.getElementById('tt');
let rate=1,playing=false;

function vids(){{return [...document.querySelectorAll('.pane[data-on] video')]}}
function attach(v,t){{
  if(!v.src&&v.dataset.src){{v.src=v.dataset.src;v.preload='auto';v.load();}}
  v.playbackRate=rate;
  // only the seek needs an event, and it needs metadata rather than data. STARTING is left
  // to the watchdog: `canplay` fires once, so a listener added after it has gone by never
  // runs -- exactly what happens when a pane opens (fires) and play is pressed a moment later.
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
  // start anything that should be running and is not -- a clip still buffering when play was
  // pressed, or one whose pane opened afterwards. Self-healing, so no event can be missed.
  if(playing)vs.forEach(v=>{{if(v.paused&&v.readyState>=2)v.play().catch(()=>{{}})}});
  const v=vs[0]; if(!v||!v.duration)return;
  sk.value=Math.round(1000*v.currentTime/v.duration);
  tt.textContent=v.currentTime.toFixed(1)+'s';
  // Only nudge clips that are ALREADY RUNNING. Seeking one that has not started restarts
  // its buffering, so a heavy clip gets yanked back every 200 ms and can never catch up --
  // it sits pinned at 0 while the light ones run away. Let it start on its own first; the
  // next tick will align it.
  vs.slice(1).forEach(o=>{{if(!o.paused&&o.duration&&
    Math.abs(o.currentTime-v.currentTime)>0.25) o.currentTime=v.currentTime}});
}},200);
toggle(tb[0].dataset.g);
</script></body></html>
"""
    out = os.path.join(a.out, "index.html")
    with open(out, "w") as fh:
        fh.write(html)
    nclip = sum(len(c) for _, rows in panes for _, c in rows)
    print(f"wrote {out} ({len(html)/1024:.0f} KB)  mode={a.mode}  "
          f"{len(panes)} tau_m, {nclip} cells, {len(imgs)} analysis plots")


if __name__ == "__main__":
    main()
