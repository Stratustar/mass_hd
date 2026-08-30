#!/usr/bin/env python3
"""A filterable video dashboard for a confluent-wet scan.

Buttons on the four things the campaign varies, keyed on the RATIOS rather than on absolute
step counts, because the ratios are what the physics is a function of:

    zeta                  the activity
    tau_m / tau_motion    how many defect encounters the memory averages over
    tau_chi / tau_motion  how far the phenotype lags behind the flow
    switch-sign           the feedback branch
    Dbio                  present as a fifth axis because without it a combination is not
                          unique -- and because the Dbio = 0 arm is grid-noise dominated,
                          which the caption says outright rather than leaving to be found

tau_motion = d/u_rms is the per-activity clock measured open-loop (20260828/cw_pmem), so a
ratio is a fixed number per (activity, tau) pair rather than something re-derived per run.

RATIO BUTTONS ARE BINNED to the {1,3} x 10^n ladder, because the raw values are near-misses
of it -- 0.095, 0.286, 0.997, 3.021, 10.574 -- and a button per raw value would be eighteen
buttons that all look the same. The bin is the label; the caption carries the exact value.
Binning is safe here only because it is checked: within one activity every (tau_m, tau_chi)
pair lands in a distinct pair of bins, so a selection still names exactly one run.

Written separately from make_dashboard.py rather than retrofitted onto it: that one keys on
absolute parameter tokens parsed out of directory names, cannot parse the sign token here
('sm' has no numeric tail), expects `<field>_<a>-<b>_step<k>.mp4` where these are plain
`chi.mp4`, and carries deliberate uncommitted local edits that should not be disturbed.

Writes <resultsdir>/index.html and filelist.txt (every referenced file), so the browsable set
travels with one `rsync -a --files-from=filelist.txt`.

Usage: cw_html.py <resultsdir> [--calib cw_calib.json] [--title T]
"""
import argparse
import json
import math
import os
import re

# Two scan namings, and they need different axes rather than one regex stretched over both.
# 20260828 cw_loop:  A9_sp_D0p03_tm330_tc1   -- absolute tau_m, tc is the tau_chi/tau_m factor,
#                    and zeta and Dbio vary, so the ratio axes have to be DERIVED and binned.
# 20260830 cw_scan:  sp_p75_tm0p2_tc0p2      -- zeta and Dbio are frozen, pmem is a percentile,
#                    and the directory already carries both ratios, so nothing is derived and
#                    nothing is binned: the button value is the number in the name.
VARIANT = re.compile(r"^A(?P<A>\d+)_s(?P<s>[pm])_D(?P<D>[\dp]+)_tm(?P<tm>\d+)_tc(?P<tc>\d+)$")
VARIANT_RATIO = re.compile(r"^s(?P<s>[pm])_p(?P<pct>\d+)_tm(?P<rm>[\dp]+)_tc(?P<rc>[\dp]+)$")
AXES_RATIO = [("rm", "tau_m / tau_motion"),
              ("rc", "tau_chi / tau_motion"),
              ("s",  "switch-sign"),
              ("pct", "pmem percentile of P")]
AXES = [("zeta", "zeta  (activity)"),
        ("rm",   "tau_m / tau_motion"),
        ("rc",   "tau_chi / tau_motion"),
        ("s",    "switch-sign"),
        ("D",    "Dbio")]
CAMPAIGN = ["collapse_montage.png", "collapse_scaled.png", "grid_power.png",
            "feedback.png", "loop_lag.png", "collapse.png"]
# tau_motion = d/u_rms per activity, measured open-loop at L=400 (20260828/cw_pmem)
TAU_MOTION = {1: 3492.0, 3: 1046.0, 9: 331.0}
CC = 0.03


def ladder_bin(x):
    """Nearest rung of the {1,3} x 10^n ladder, in the geometric sense.

    The rungs are 1 and 3 per decade, so the boundaries sit at their geometric means,
    sqrt(3) = 1.73 and sqrt(30) = 5.48 -- above 5.48 a value belongs to the NEXT decade's 1,
    not to this decade's 3. Missing that carries 0.095 to 0.03 instead of 0.1 and 0.859 to
    0.3 instead of 1, which then collides two distinct runs into one button pair.
    """
    if x <= 0:
        return 0.0
    e = math.floor(math.log10(x))
    m = x / 10.0**e
    if m < math.sqrt(3.0):
        return 10.0**e
    if m < math.sqrt(30.0):
        return 3.0 * 10.0**e
    return 10.0**(e + 1)


def fmt(x):
    return f"{x:g}"


def load_tau_motion(path):
    if not path or not os.path.exists(path):
        return dict(TAU_MOTION)
    with open(path) as f:
        rows = json.load(f)
    out = {}
    for r in rows:
        try:
            out[int(round(r["A"]))] = float(r["tau_motion"])
        except Exception:
            pass
    return out or dict(TAU_MOTION)


def parse(name, tmot):
    mt = VARIANT_RATIO.match(name)
    if mt:
        return {"rm": fmt(float(mt["rm"].replace("p", "."))),
                "rc": fmt(float(mt["rc"].replace("p", "."))),
                "s": "+1" if mt["s"] == "p" else "-1",
                "pct": mt["pct"],
                "_A": None, "_tm": None, "_tc": None,
                "_rm": float(mt["rm"].replace("p", ".")),
                "_rc": float(mt["rc"].replace("p", ".")), "_tmot": tmot.get(9)}
    mt = VARIANT.match(name)
    if not mt:
        return None
    A = int(mt["A"])
    tm, k = float(mt["tm"]), float(mt["tc"])
    tc = tm * k
    t0 = tmot.get(A)
    if not t0:
        return None
    rm, rc = tm / t0, tc / t0
    return {"zeta": fmt(round(2 * A * CC, 4)),          # zeta = zeta_eff/(1-chi0), chi0 = .5
            "rm": fmt(ladder_bin(rm)), "rc": fmt(ladder_bin(rc)),
            "s": "+1" if mt["s"] == "p" else "-1",
            "D": mt["D"].replace("p", "."),
            "_A": A, "_tm": tm, "_tc": tc, "_rm": rm, "_rc": rc, "_tmot": t0}


# THE standard order, matching cw_dash.py. Changing it changes the standard.
VIDEOS = ["u", "pressure", "m", "chi"]
VIDEO_DIR = "dash"
ALL_FIGS = {"crosscorr": "C(r), cross",
            "autocorr": "C(r), auto",
            "lagcorr": "lagged cross-correlation",
            "st": "chi spatio-temporal C(r,tau)"}
DEFAULT_FIGS = ["crosscorr", "st"]


def collect(root, tmot, figs):
    """The four clean field videos and the three correlation figures, and nothing else.

    Everything the diagnostic build carried -- the dashboard video, the six-panel sheet, the
    scalar series, the field stills -- is still on disk and is deliberately not listed. A page
    that shows everything is a page nobody reads.
    """
    cases = {}
    for d in sorted(os.listdir(root)):
        if not os.path.isdir(os.path.join(root, d)):
            continue
        p = parse(d, tmot)
        if p is None:
            continue
        p["_v"] = [(n, f"{d}/{VIDEO_DIR}/{n}.mp4") for n in VIDEOS
                   if os.path.exists(os.path.join(root, d, VIDEO_DIR, f"{n}.mp4"))]
        p["_f"] = [(ALL_FIGS.get(t, t), f"clean/{d}_{t}.png") for t in figs
                   if os.path.exists(os.path.join(root, "clean", f"{d}_{t}.png"))]
        cases[d] = p
    return cases


def captions(root, cases):
    """Only the run's name. The buttons already say what the parameters are."""
    return {n: n for n in cases}


HTML = """<!doctype html><meta charset="utf-8"><title>%(title)s</title>
<style>
 body{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;margin:0;
      background:#fff;color:#18181b}
 .filters{padding:14px 26px 16px;border-bottom:1px solid #ececef;position:sticky;top:0;
          background:#fff;z-index:5}
 .row{display:flex;align-items:center;gap:7px;margin:6px 0;flex-wrap:wrap}
 .row b{min-width:170px;font-weight:600;color:#3f3f46;font-size:13px}
 button{border:1px solid #d4d4d8;background:#fff;border-radius:6px;padding:5px 13px;
        cursor:pointer;font-size:13px;font-variant-numeric:tabular-nums}
 button.on{background:#1d4ed8;border-color:#1d4ed8;color:#fff}
 button:disabled{opacity:.26;cursor:not-allowed}
 main{padding:20px 26px 40px}
 .name{color:#a1a1aa;font:12px ui-monospace,SFMono-Regular,monospace;margin:0 0 14px}
 .vids{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-bottom:26px}
 .figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(460px,1fr));gap:18px}
 /* The four field videos are a ROW, not a grid: the standard is u, p, m, chi read left to
    right, and wrapping them to 2x2 at a 1400px breakpoint broke that on any normal laptop
    pane. They stay four-across and simply get smaller; only a genuinely narrow window
    (<900px, where each would be under 210px) is allowed to wrap. */
 @media(max-width:900px){.vids{grid-template-columns:repeat(2,1fr)}}
 @media(max-width:560px){.vids{grid-template-columns:1fr}}
 figure{margin:0}
 video,img{width:100%%;height:auto;display:block;border-radius:5px}
 .none{color:#a1a1aa;padding:34px 0}
</style>
<div class="filters" id="f"></div>
<main><p class="name" id="cap"></p>
<div class="vids" id="v"></div><div class="figs" id="g"></div></main>
<script>
const CASES=%(cases)s, CAPS=%(caps)s, AXES=%(axes)s;
const ORDER={}; AXES.forEach(([k])=>{
  ORDER[k]=[...new Set(Object.values(CASES).map(c=>c[k]))].sort((a,b)=>parseFloat(a)-parseFloat(b));});
const sel={}; AXES.forEach(([k])=>{sel[k]=ORDER[k][Math.min(2,ORDER[k].length-1)];});
const match=(c,skip)=>AXES.every(([k])=>k===skip||c[k]===sel[k]);
function draw(){
  const f=document.getElementById('f'); f.innerHTML='';
  AXES.forEach(([k,label])=>{
    const row=document.createElement('div'); row.className='row';
    row.innerHTML='<b>'+label+'</b>';
    ORDER[k].forEach(v=>{
      const b=document.createElement('button'); b.textContent=v;
      if(sel[k]===v) b.className='on';
      else if(!Object.values(CASES).some(c=>c[k]===v&&match(c,k))) b.className='dim';
      // A dimmed value is still clickable. The two tau axes are coupled -- tau_chi >= tau_m
      // always -- so reaching some runs means moving both at once, and a disabled button
      // makes those runs unreachable: you cannot take the first step of a two-step move.
      // Clicking one instead snaps the OTHER axes to the nearest run that has this value.
      b.onclick=()=>{
        const cand=Object.entries(CASES).filter(([,c])=>c[k]===v);
        if(!cand.length) return;
        let best=null, bestd=1e9;
        for(const [n,c] of cand){
          let d=0;
          for(const [j] of AXES){ if(j===k||c[j]===sel[j]) continue;
            const o=ORDER[j], a=o.indexOf(c[j]), b2=o.indexOf(sel[j]);
            d += 1 + Math.abs(a-b2)/o.length;   // prefer the smallest move, on fewest axes
          }
          if(d<bestd){bestd=d; best=c;}
        }
        AXES.forEach(([j])=>{sel[j]=best[j];});
        sel[k]=v;
        draw();
      };
      row.appendChild(b);
    });
    f.appendChild(row);
  });
  const name=Object.keys(CASES).find(n=>match(CASES[n]));
  const v=document.getElementById('v'), g=document.getElementById('g');
  v.innerHTML=''; g.innerHTML='';
  document.getElementById('cap').textContent = name ? (CAPS[name]||name) : '';
  if(!name){ v.innerHTML='<div class="none">no run at this combination</div>'; return; }
  const c=CASES[name];
  (c._v||[]).forEach(([n,src])=>{
    const d=document.createElement('figure');
    d.innerHTML='<video src="'+src+'" controls autoplay loop muted playsinline></video>';
    v.appendChild(d);
  });
  (c._f||[]).forEach(([lab,src])=>{
    const d=document.createElement('figure');
    d.innerHTML='<img loading="lazy" alt="'+lab+'" src="'+src+'">';
    g.appendChild(d);
  });
}
draw();
</script>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resultsdir")
    ap.add_argument("--calib", default=None,
                    help="cw_calib.json supplying tau_motion per activity")
    ap.add_argument("--title", default="confluent-wet_closed-loop_scan")
    ap.add_argument("--figs", nargs="*", default=None,
                    help=f"static figures to show, in order. default {DEFAULT_FIGS}; "
                         f"available {sorted(ALL_FIGS)}")
    a = ap.parse_args()
    root = a.resultsdir

    tmot = load_tau_motion(a.calib)
    cases = collect(root, tmot, a.figs or DEFAULT_FIGS)
    if cases and "rm" in next(iter(cases.values())) and next(iter(cases.values()))["_A"] is None:
        global AXES
        AXES = AXES_RATIO
    if not cases:
        raise SystemExit(f"no variant directories under {root}")

    # the binning is only legitimate if it still names a unique run -- check, do not assume
    keys = {}
    for n, c in cases.items():
        k = tuple(c[ax] for ax, _ in AXES)
        keys.setdefault(k, []).append(n)
    clashes = {k: v for k, v in keys.items() if len(v) > 1}
    if clashes:
        print(f"!! {len(clashes)} button combinations map to more than one run:")
        for k, v in list(clashes.items())[:5]:
            print("   ", k, "->", v)
    else:
        print(f"binning check: all {len(cases)} runs are uniquely addressable by the buttons")

    caps = captions(root, cases)
    referenced = []
    for c in cases.values():
        referenced += [p for _, p in c["_v"]] + [p for _, p in c["_f"]]

    slim = {n: {**{k: c[k] for k, _ in AXES}, "_v": c["_v"], "_f": c["_f"]}
            for n, c in cases.items()}
    with open(os.path.join(root, "index.html"), "w") as f:
        f.write(HTML % dict(title=a.title, cases=json.dumps(slim),
                            caps=json.dumps(caps), axes=json.dumps(AXES)))
    referenced = sorted(set(p for p in referenced if os.path.exists(os.path.join(root, p))))
    with open(os.path.join(root, "filelist.txt"), "w") as f:
        f.write("\n".join(["index.html"] + referenced) + "\n")
    mb = sum(os.path.getsize(os.path.join(root, p)) for p in referenced) / 1e6
    print(f"wrote {os.path.join(root, 'index.html')}")
    nv = sum(len(c["_v"]) for c in cases.values())
    nf = sum(len(c["_f"]) for c in cases.values())
    print(f"  {len(cases)} runs, {nv} videos, {nf} figures, {mb:.0f} MB")
    for ax, lab in AXES:
        vals = sorted({c[ax] for c in cases.values()}, key=float)
        print(f"  {lab:24s}: {' '.join(vals)}")
    print(f"  rsync -a --files-from=filelist.txt jed:{root}/ <local>/")


if __name__ == "__main__":
    main()
