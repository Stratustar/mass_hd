#!/usr/bin/env python3
"""A self-contained HTML index for browsing a confluent-wet scan: 72 cases, 5 parameters.

Same idea as make_dashboard.py -- per-parameter filter buttons with unavailable combinations
greyed out -- but written separately rather than retrofitted onto that one, because the
confluent-wet outputs differ in every way it keys on: the variant names carry a sign token
(sp/sm) that its TOKEN_RE cannot parse, the videos are plain `<field>.mp4` instead of
`<field>_<a>-<b>_step<k>.mp4`, and make_dashboard.py currently carries deliberate
uncommitted local edits that should not be disturbed.

Writes <resultsdir>/index.html and <resultsdir>/filelist.txt, the latter being every file the
page references, so the whole thing can be pulled down with

    rsync -a --files-from=filelist.txt jed:<resultsdir>/ <local>/

Usage: cw_html.py <resultsdir> [--title T]
"""
import argparse
import json
import os
import re

VARIANT = re.compile(r"^A(?P<A>\d+)_s(?P<s>[pm])_D(?P<D>[\dp]+)_tm(?P<tm>\d+)_tc(?P<tc>\d+)$")
AXES = [("A", "activity A"), ("s", "switch-sign"), ("D", "Dbio"),
        ("tm", "tau_m"), ("tc", "tau_chi / tau_m")]
CAMPAIGN = ["collapse_montage.png", "collapse_scaled.png", "grid_power.png",
            "feedback.png", "loop_lag.png", "collapse.png"]


def num(v):
    return v.replace("p", ".")


def parse(name):
    mt = VARIANT.match(name)
    if not mt:
        return None
    return {"A": mt["A"], "s": ("+1" if mt["s"] == "p" else "-1"),
            "D": num(mt["D"]), "tm": mt["tm"], "tc": mt["tc"]}


def collect(root):
    cases = {}
    for d in sorted(os.listdir(root)):
        p = parse(d)
        if p is None or not os.path.isdir(os.path.join(root, d)):
            continue
        files = sorted(os.listdir(os.path.join(root, d)))
        media = {
            "dashboard": [f for f in files if f == "dashboard.mp4"],
            "chi": [f for f in files if f == "chi.mp4"],
            "stills": [f for f in files if f.startswith("fields_") and f.endswith(".png")],
            "scalars": [f for f in files if f == "scalars_vs_t.png"],
        }
        sheet = os.path.join("figs", f"case_{d}.png")
        p["_files"] = media
        p["_sheet"] = sheet if os.path.exists(os.path.join(root, sheet)) else None
        cases[d] = p
    return cases


def metrics(root):
    """A few headline numbers per case, if the merged scan is there, for the caption line."""
    path = os.path.join(root, "cw_scan.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        rows = json.load(f)
    out = {}
    for r in rows:
        try:
            out[r["case"]] = (f"A_eff = {r['A_eff']:.2f}   L_chi/d = {r['L_chi']/r['d']:.2f}   "
                              f"tau_m/t_eddy = {r['tau_m_set']/r['t_eddy']:.2f}   "
                              f"grid-power(chi) = {r['grid_power_chi']:.1e}")
        except Exception:
            pass
    return out


HTML = """<!doctype html><meta charset="utf-8"><title>%(title)s</title>
<style>
 body{font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;margin:0;
      background:#fafafa;color:#1a1a1a}
 header{padding:18px 24px;background:#fff;border-bottom:1px solid #e3e3e3}
 h1{margin:0 0 4px;font-size:19px} .sub{color:#666;font-size:13px}
 .filters{padding:14px 24px;background:#fff;border-bottom:1px solid #e3e3e3}
 .row{display:flex;align-items:center;gap:8px;margin:6px 0;flex-wrap:wrap}
 .row b{min-width:120px;font-weight:600;color:#444;font-size:13px}
 button{border:1px solid #ccc;background:#fff;border-radius:5px;padding:4px 11px;
        cursor:pointer;font-size:13px}
 button.on{background:#1a5fb4;border-color:#1a5fb4;color:#fff}
 button:disabled{opacity:.3;cursor:not-allowed}
 main{padding:20px 24px} .cap{color:#555;margin:0 0 10px;font-family:ui-monospace,monospace}
 .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:16px}
 figure{margin:0;background:#fff;border:1px solid #e3e3e3;border-radius:7px;padding:10px}
 figcaption{font-size:12px;color:#666;margin-bottom:6px}
 video,img{width:100%%;height:auto;border-radius:4px;display:block}
 .wide{grid-column:1/-1}
 details{margin-top:22px;background:#fff;border:1px solid #e3e3e3;border-radius:7px;padding:12px}
 summary{cursor:pointer;font-weight:600}
 .none{color:#999;padding:30px 0}
</style>
<header><h1>%(title)s</h1>
<div class="sub">%(ncase)d cases &middot; filter buttons grey out combinations that were not run</div></header>
<div class="filters" id="f"></div>
<main><div class="cap" id="cap"></div><div class="grid" id="g"></div>
<details><summary>campaign figures</summary><div class="grid">%(campaign)s</div></details></main>
<script>
const CASES = %(cases)s, METRICS = %(metrics)s;
const AXES = %(axes)s;
const sel = {};
AXES.forEach(([k])=>{ sel[k] = [...new Set(Object.values(CASES).map(c=>c[k]))]
  .sort((a,b)=>parseFloat(a)-parseFloat(b))[0]; });
function match(c, skip){ return AXES.every(([k])=> k===skip || c[k]===sel[k]); }
function draw(){
  const f = document.getElementById('f'); f.innerHTML='';
  AXES.forEach(([k,label])=>{
    const row=document.createElement('div'); row.className='row';
    row.innerHTML='<b>'+label+'</b>';
    const vals=[...new Set(Object.values(CASES).map(c=>c[k]))]
      .sort((a,b)=>parseFloat(a)-parseFloat(b));
    vals.forEach(v=>{
      const b=document.createElement('button'); b.textContent=v;
      if(sel[k]===v) b.className='on';
      const ok=Object.values(CASES).some(c=>c[k]===v && match(c,k));
      b.disabled=!ok && sel[k]!==v;
      b.onclick=()=>{sel[k]=v;draw();};
      row.appendChild(b);
    });
    f.appendChild(row);
  });
  const name = Object.keys(CASES).find(n=>match(CASES[n]));
  const g=document.getElementById('g'); g.innerHTML='';
  document.getElementById('cap').textContent = name ? (name+'    '+(METRICS[name]||'')) : '';
  if(!name){ g.innerHTML='<div class="none">no run at this combination</div>'; return; }
  const c=CASES[name], add=(html,wide)=>{
    const d=document.createElement('figure'); if(wide) d.className='wide';
    d.innerHTML=html; g.appendChild(d); };
  (c._files.dashboard||[]).forEach(f=>add('<figcaption>dashboard &mdash; chi, m, P, |u| with the '
    +'loop time series</figcaption><video src="'+name+'/'+f+'" controls autoplay loop muted></video>',1));
  (c._files.chi||[]).forEach(f=>add('<figcaption>chi</figcaption><video src="'+name+'/'+f
    +'" controls loop muted></video>'));
  if(c._sheet) add('<figcaption>correlations: spatial, lagged, C(r,tau), series, spectra'
    +'</figcaption><img loading="lazy" src="'+c._sheet+'">',1);
  (c._files.scalars||[]).forEach(f=>add('<figcaption>scalars vs t</figcaption>'
    +'<img loading="lazy" src="'+name+'/'+f+'">'));
  (c._files.stills||[]).forEach(f=>add('<figcaption>'+f+'</figcaption>'
    +'<img loading="lazy" src="'+name+'/'+f+'">'));
}
draw();
</script>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resultsdir")
    ap.add_argument("--title", default="confluent-wet closed-loop scan")
    a = ap.parse_args()
    root = a.resultsdir

    cases = collect(root)
    if not cases:
        raise SystemExit(f"no variant directories under {root}")
    met = metrics(root)

    referenced = []
    for name, c in cases.items():
        for grp in c["_files"].values():
            referenced += [os.path.join(name, f) for f in grp]
        if c["_sheet"]:
            referenced.append(c["_sheet"])
    camp_html, camp_files = "", []
    for f in CAMPAIGN:
        p = os.path.join("figs", f)
        if os.path.exists(os.path.join(root, p)):
            camp_html += (f'<figure class="wide"><figcaption>{f}</figcaption>'
                          f'<img loading="lazy" src="{p}"></figure>')
            camp_files.append(p)
    referenced += camp_files

    slim = {n: {k: c[k] for k, _ in AXES} | {"_files": c["_files"], "_sheet": c["_sheet"]}
            for n, c in cases.items()}
    html = HTML % dict(title=a.title, ncase=len(cases), campaign=camp_html,
                       cases=json.dumps(slim), metrics=json.dumps(met),
                       axes=json.dumps(AXES))
    out = os.path.join(root, "index.html")
    with open(out, "w") as f:
        f.write(html)
    with open(os.path.join(root, "filelist.txt"), "w") as f:
        f.write("\n".join(["index.html"] + sorted(set(referenced))) + "\n")
    mb = sum(os.path.getsize(os.path.join(root, p)) for p in set(referenced)
             if os.path.exists(os.path.join(root, p))) / 1e6
    print(f"wrote {out}\n  {len(cases)} cases, {len(set(referenced))} files, {mb:.0f} MB")
    print(f"  pull with:  rsync -a --files-from=filelist.txt jed:{root}/ <local>/")


if __name__ == "__main__":
    main()
