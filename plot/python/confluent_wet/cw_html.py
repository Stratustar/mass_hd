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

VARIANT = re.compile(r"^A(?P<A>\d+)_s(?P<s>[pm])_D(?P<D>[\dp]+)_tm(?P<tm>\d+)_tc(?P<tc>\d+)$")
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


def collect(root, tmot):
    cases = {}
    for d in sorted(os.listdir(root)):
        if not os.path.isdir(os.path.join(root, d)):
            continue
        p = parse(d, tmot)
        if p is None:
            continue
        files = sorted(os.listdir(os.path.join(root, d)))
        p["_v"] = [f for f in ("dashboard.mp4", "chi.mp4") if f in files]
        p["_img"] = ([os.path.join("figs", f"case_{d}.png")]
                     if os.path.exists(os.path.join(root, "figs", f"case_{d}.png")) else [])
        p["_img"] += [f for f in files if f == "scalars_vs_t.png"]
        p["_img"] += [f for f in files if f.startswith("fields_") and f.endswith(".png")]
        cases[d] = p
    return cases


def captions(root, cases):
    path = os.path.join(root, "cw_scan.json")
    scan = {}
    if os.path.exists(path):
        with open(path) as f:
            for r in json.load(f):
                scan[r["case"]] = r
    out = {}
    for n, c in cases.items():
        base = (f"A = {c['_A']}   zeta = {c['zeta']}   tau_motion = {c['_tmot']:.0f}   "
                f"tau_m = {c['_tm']:.0f} ({c['_rm']:.2f} x tau_motion)   "
                f"tau_chi = {c['_tc']:.0f} ({c['_rc']:.2f} x tau_motion)")
        r = scan.get(n)
        if r:
            base += (f"\nA_eff = {r['A_eff']:.2f}   L_chi/d = {r['L_chi']/r['d']:.2f}   "
                     f"L_chi = {r['L_chi']:.1f}   d = {r['d']:.1f}   "
                     f"<g(P)> = {r['g_bar']:.2f}   grid-power(chi) = {r['grid_power_chi']:.1e}")
            if r["grid_power_chi"] > 1e-3:
                base += "   <-- GRID-NOISE DOMINATED, structure not usable"
        out[n] = base
    return out


HTML = """<!doctype html><meta charset="utf-8"><title>%(title)s</title>
<style>
 body{font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;margin:0;
      background:#f7f7f8;color:#18181b}
 header{padding:16px 24px;background:#fff;border-bottom:1px solid #e4e4e7}
 h1{margin:0 0 3px;font-size:19px} .sub{color:#71717a;font-size:13px}
 .filters{padding:12px 24px 16px;background:#fff;border-bottom:1px solid #e4e4e7;
          position:sticky;top:0;z-index:5}
 .row{display:flex;align-items:center;gap:7px;margin:7px 0;flex-wrap:wrap}
 .row b{min-width:175px;font-weight:600;color:#3f3f46;font-size:13px}
 button{border:1px solid #d4d4d8;background:#fff;border-radius:6px;padding:5px 13px;
        cursor:pointer;font-size:13px;font-variant-numeric:tabular-nums}
 button.on{background:#1d4ed8;border-color:#1d4ed8;color:#fff}
 button:disabled{opacity:.28;cursor:not-allowed}
 main{padding:18px 24px}
 .cap{color:#3f3f46;margin:0 0 12px;font:12.5px/1.6 ui-monospace,SFMono-Regular,monospace;
      white-space:pre-wrap;background:#fff;border:1px solid #e4e4e7;border-radius:7px;
      padding:10px 13px}
 .warn{color:#b91c1c;font-weight:600}
 .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:15px}
 figure{margin:0;background:#fff;border:1px solid #e4e4e7;border-radius:8px;padding:11px}
 figcaption{font-size:12px;color:#71717a;margin-bottom:7px}
 video,img{width:100%%;height:auto;border-radius:5px;display:block}
 .wide{grid-column:1/-1}
 details{margin-top:22px;background:#fff;border:1px solid #e4e4e7;border-radius:8px;padding:13px}
 summary{cursor:pointer;font-weight:600}
 .none{color:#a1a1aa;padding:34px 0}
</style>
<header><h1>%(title)s</h1>
<div class="sub">%(ncase)d runs &middot; ratio buttons are binned to the {1,3}&times;10<sup>n</sup>
 ladder; the exact values are in the caption &middot; combinations that were not run are greyed out</div></header>
<div class="filters" id="f"></div>
<main><div class="cap" id="cap"></div><div class="grid" id="g"></div>
<details><summary>campaign figures</summary><div class="grid">%(campaign)s</div></details></main>
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
      b.disabled=!Object.values(CASES).some(c=>c[k]===v&&match(c,k))&&sel[k]!==v;
      b.onclick=()=>{
        sel[k]=v;
        // keep the selection on a run that exists: relax the other axes if needed
        if(!Object.keys(CASES).some(n=>match(CASES[n]))){
          for(const [j] of AXES){ if(j===k) continue;
            const alt=ORDER[j].find(w=>{const s=sel[j];sel[j]=w;
              const ok=Object.keys(CASES).some(n=>match(CASES[n]));sel[j]=s;return ok;});
            if(alt!==undefined){sel[j]=alt; if(Object.keys(CASES).some(n=>match(CASES[n])))break;}
          }
        }
        draw();
      };
      row.appendChild(b);
    });
    f.appendChild(row);
  });
  const name=Object.keys(CASES).find(n=>match(CASES[n]));
  const g=document.getElementById('g'); g.innerHTML='';
  const cap=document.getElementById('cap');
  if(!name){ cap.textContent=''; g.innerHTML='<div class="none">no run at this combination</div>'; return; }
  const txt=name+'\\n'+(CAPS[name]||'');
  cap.innerHTML=txt.replace(/(<-- GRID-NOISE DOMINATED[^\\n]*)/,'<span class="warn">$1</span>');
  const c=CASES[name];
  const add=(h,w)=>{const d=document.createElement('figure');if(w)d.className='wide';
                    d.innerHTML=h;g.appendChild(d);};
  (c._v||[]).forEach(v=>{
    const lab = v==='dashboard.mp4'
      ? 'dashboard &mdash; chi, m, P, |u| with the loop time series'
      : 'chi';
    add('<figcaption>'+lab+'</figcaption><video src="'+name+'/'+v+'" controls loop muted '
        +(v==='dashboard.mp4'?'autoplay ':'')+'></video>', v==='dashboard.mp4');
  });
  (c._img||[]).forEach(p=>{
    const src = p.indexOf('/')>=0 ? p : name+'/'+p;
    add('<figcaption>'+p.split('/').pop()+'</figcaption><img loading="lazy" src="'+src+'">',
        p.indexOf('figs/')===0);
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
    a = ap.parse_args()
    root = a.resultsdir

    tmot = load_tau_motion(a.calib)
    cases = collect(root, tmot)
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
    for n, c in cases.items():
        referenced += [os.path.join(n, f) for f in c["_v"]]
        referenced += [p if p.startswith("figs/") else os.path.join(n, p) for p in c["_img"]]
    camp_html = ""
    for f in CAMPAIGN:
        p = os.path.join("figs", f)
        if os.path.exists(os.path.join(root, p)):
            camp_html += (f'<figure class="wide"><figcaption>{f}</figcaption>'
                          f'<img loading="lazy" src="{p}"></figure>')
            referenced.append(p)

    slim = {n: {**{k: c[k] for k, _ in AXES}, "_v": c["_v"], "_img": c["_img"]}
            for n, c in cases.items()}
    with open(os.path.join(root, "index.html"), "w") as f:
        f.write(HTML % dict(title=a.title, ncase=len(cases), campaign=camp_html,
                            cases=json.dumps(slim), caps=json.dumps(caps),
                            axes=json.dumps(AXES)))
    referenced = sorted(set(p for p in referenced if os.path.exists(os.path.join(root, p))))
    with open(os.path.join(root, "filelist.txt"), "w") as f:
        f.write("\n".join(["index.html"] + referenced) + "\n")
    mb = sum(os.path.getsize(os.path.join(root, p)) for p in referenced) / 1e6
    print(f"wrote {os.path.join(root, 'index.html')}")
    print(f"  {len(cases)} runs, {len(referenced)} files, {mb:.0f} MB")
    for ax, lab in AXES:
        vals = sorted({c[ax] for c in cases.values()}, key=float)
        print(f"  {lab:24s}: {' '.join(vals)}")
    print(f"  rsync -a --files-from=filelist.txt jed:{root}/ <local>/")


if __name__ == "__main__":
    main()
