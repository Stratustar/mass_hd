#!/usr/bin/env python3
"""Generate a self-contained HTML dashboard to browse simulation result figures.

Scans a study results directory whose subdirectories are parameter-encoded
variants (e.g. a0p0001_z0p05_O0p005_LL0p032_xi1p2), builds clickable
per-parameter filter buttons, and shows the gifs + last-frame static figures
for the selected parameter combination. Click any figure/gif to enlarge.

Usage:
    python make_dashboard.py <study_results_dir> [--title T] [--out index.html]

Writes <study_results_dir>/<out> and <study_results_dir>/dashboard_filelist.txt
(relative paths of every referenced file + the html, for `rsync --files-from`).
"""
import argparse
import json
import os
import re

PARAM_ORDER = ["alpha", "zeta", "Ochi", "LL", "xi", "Dchi", "Achi", "friction"]
KEY_LABEL = {"a": "alpha", "z": "zeta", "O": "Ochi", "LL": "LL", "xi": "xi",
             "Dchi": "Dchi", "Achi": "Achi", "fr": "friction"}

# Figure groups not shown in the dashboard (files stay on disk, just not listed).
EXCLUDE_GROUPS = {
    "field_snapshots", "final_state_diagnostics", "mask_diagnostics",
    "proliferation_dashboard", "radial_kymographs", "radial_profiles",
}

GIF_RE = re.compile(r"^(?P<field>.+)_(?P<a>\d+)-(?P<b>\d+)_step(?P<k>\d+)\.gif$")
FRAME_RE = re.compile(r"^(?P<field>.+)_frame(?P<n>\d+)\.png$")
TOKEN_RE = re.compile(r"^([A-Za-z]+)([-0-9pe.]+)$")


def parse_variant(name):
    params = {}
    for tok in name.split("_"):
        m = TOKEN_RE.match(tok)
        if not m:
            continue
        key, raw = m.group(1), m.group(2)
        params[KEY_LABEL.get(key, key)] = raw.replace("p", ".")
    # Show Achi as a multiple of Ochi (scans set Achi = k*Ochi), so the filter
    # button reads "0.5xOchi" / "1xOchi" instead of the absolute value.
    if "Achi" in params and "Ochi" in params:
        try:
            a, o = float(params["Achi"]), float(params["Ochi"])
            if o != 0:
                params["Achi"] = "%gxOchi" % (a / o)
        except (TypeError, ValueError):
            pass
    return params


def collect_figures(vdir):
    try:
        files = os.listdir(vdir)
    except OSError:
        return []
    gifs, frames, singles = [], {}, []
    for f in files:
        gm = GIF_RE.match(f)
        if gm:
            gifs.append((gm.group("field"), f))
            continue
        fm = FRAME_RE.match(f)
        if fm:
            field, n = fm.group("field"), int(fm.group("n"))
            if field not in frames or n > frames[field][0]:
                frames[field] = (n, f)
            continue
        if f.endswith(".png"):
            singles.append(f)
    last_png = {field: f for field, (n, f) in frames.items()}
    single_groups = {f[:-4] for f in singles}
    figs = []
    for field, f in sorted(gifs):
        fig = {"kind": "gif", "group": field, "label": field + " (gif)", "file": f}
        if field in last_png:
            fig["poster"] = last_png[field]   # lightweight cover; gif itself loads on click
        figs.append(fig)
    for field, (n, f) in sorted(frames.items()):
        # if a same-named single (e.g. a pressure decomposition pressure.png) exists, it
        # replaces the plain last-frame tile; the frame is kept only as the gif poster.
        if field in single_groups:
            continue
        figs.append({"kind": "png", "group": field, "label": "%s (frame %d)" % (field, n), "file": f})
    for f in sorted(singles):
        figs.append({"kind": "png", "group": f[:-4], "label": f[:-4], "file": f})
    return [f for f in figs if f["group"] not in EXCLUDE_GROUPS]


def sort_key(s):
    try:
        return (0, float(s))
    except (TypeError, ValueError):
        return (1, str(s))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("studydir")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default="index.html")
    args = ap.parse_args()

    study = os.path.abspath(args.studydir)
    title = args.title or os.path.basename(study)

    variants, referenced = [], []
    for name in sorted(os.listdir(study)):
        vdir = os.path.join(study, name)
        if not os.path.isdir(vdir):
            continue
        figs = collect_figures(vdir)
        if not figs:
            continue
        variants.append({"name": name, "params": parse_variant(name), "figures": figs})
        for fig in figs:
            referenced.append("%s/%s" % (name, fig["file"]))

    axes = {}
    for v in variants:
        for k, val in v["params"].items():
            axes.setdefault(k, set()).add(val)
    axes = {k: sorted(vals, key=sort_key) for k, vals in axes.items()}
    order = [k for k in PARAM_ORDER if k in axes] + [k for k in sorted(axes) if k not in PARAM_ORDER]

    data = {"title": title, "order": order,
            "axes": {k: axes[k] for k in order}, "variants": variants}

    out_path = os.path.join(study, args.out)
    with open(out_path, "w") as fh:
        fh.write(TEMPLATE.replace("/*DATA*/", json.dumps(data)))
    with open(os.path.join(study, "dashboard_filelist.txt"), "w") as fh:
        fh.write(args.out + "\n")
        fh.write("\n".join(referenced) + ("\n" if referenced else ""))

    print("variants=%d  figures=%d  wrote %s" % (len(variants), len(referenced), out_path))


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>results dashboard</title>
<style>
  :root { --bg:#fafafa; --fg:#1a1a1a; --muted:#666; --line:#ddd; --accent:#c0653a; --chip:#fff; }
  * { box-sizing: border-box; }
  body { margin:0; font-family: system-ui, -apple-system, sans-serif; color:var(--fg); background:var(--bg); }
  header { position:sticky; top:0; z-index:5; background:var(--bg); border-bottom:1px solid var(--line); padding:12px 16px; }
  h1 { font-size:16px; margin:0 0 8px; }
  .axis { display:flex; align-items:center; gap:6px; flex-wrap:wrap; margin:3px 0; }
  .axis .name { width:64px; font-size:12px; color:var(--muted); text-align:right; font-weight:600; }
  button.v { border:1px solid var(--line); background:var(--chip); color:var(--fg); border-radius:14px;
             padding:3px 10px; font-size:12px; cursor:pointer; font-variant-numeric: tabular-nums; }
  button.v:hover { border-color:var(--accent); }
  button.v.on { background:var(--accent); color:#fff; border-color:var(--accent); }
  button.v:disabled { opacity:.3; cursor:not-allowed; background:var(--chip); color:var(--muted);
                      border-style:dashed; border-color:var(--line); }
  button.v:disabled:hover { border-color:var(--line); }
  .bar { display:flex; align-items:center; gap:12px; margin-top:8px; flex-wrap:wrap; }
  .bar .count { font-size:12px; color:var(--muted); }
  button.act { border:1px solid var(--line); background:var(--chip); border-radius:6px; padding:3px 10px; font-size:12px; cursor:pointer; }
  .figtypes { margin-top:6px; }
  .figtypes summary { font-size:12px; color:var(--muted); cursor:pointer; }
  main { padding:14px 16px 60px; }
  .card { border:1px solid var(--line); border-radius:10px; padding:10px 12px; margin-bottom:14px; background:#fff; }
  .card h2 { font-size:13px; margin:0 0 8px; font-variant-numeric: tabular-nums; }
  .card h2 .k { color:var(--muted); font-weight:500; }
  .grid { display:grid; grid-template-columns: repeat(auto-fill, minmax(220px,1fr)); gap:10px; }
  figure { position:relative; margin:0; border:1px solid var(--line); border-radius:8px; overflow:hidden; background:#f3f3f3; }
  figure img { width:100%; height:170px; object-fit:contain; display:block; background:#7e7e7e; cursor:zoom-in; }
  .badge { position:absolute; top:6px; left:6px; background:rgba(0,0,0,.6); color:#fff; font-size:11px; padding:1px 6px; border-radius:10px; pointer-events:none; }
  .noposter { height:170px; display:flex; align-items:center; justify-content:center; color:#fff; background:#555; cursor:pointer; font-size:13px; }
  figcaption { font-size:11px; color:var(--muted); padding:4px 6px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
  .empty { color:var(--muted); font-size:13px; padding:20px 0; }
  #lb { position:fixed; inset:0; background:rgba(0,0,0,.85); display:none; align-items:center; justify-content:center; z-index:20; cursor:zoom-out; }
  #lb img { max-width:96vw; max-height:92vh; object-fit:contain; }
  #lb .cap { position:fixed; bottom:10px; left:0; right:0; text-align:center; color:#eee; font-size:13px; }
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <div id="filters"></div>
  <div class="bar">
    <button class="act" id="reset">Reset</button>
    <span class="count" id="count"></span>
  </div>
  <div id="figfilters"></div>
</header>
<main><div id="results"></div></main>
<div id="lb"><img alt=""><div class="cap"></div></div>
<script>
const DATA = /*DATA*/;
const sel = {};                 // axis -> Set of selected value strings
const axisBtns = {};            // axis -> [{val, btn}] for availability greying
const figOff = new Set();       // field groups toggled OFF (empty = all shown)
const kindOff = new Set();      // kinds toggled OFF: "static" and/or "gif"

document.getElementById("title").textContent = DATA.title;

function enc(p){ return p.split("/").map(encodeURIComponent).join("/"); }

function buildFilters(){
  const host = document.getElementById("filters");
  DATA.order.forEach(k => {
    sel[k] = new Set();
    axisBtns[k] = [];
    const row = document.createElement("div"); row.className = "axis";
    const nm = document.createElement("span"); nm.className = "name"; nm.textContent = k; row.appendChild(nm);
    DATA.axes[k].forEach(val => {
      const b = document.createElement("button"); b.className = "v"; b.textContent = val;
      b.onclick = () => { b.classList.toggle("on"); if(sel[k].has(val)) sel[k].delete(val); else sel[k].add(val); render(); };
      row.appendChild(b);
      axisBtns[k].push({val, btn: b});
    });
    host.appendChild(row);
  });
  // figure filters (visible buttons): which field groups + which kinds to show
  const groups = [];
  const seen = new Set();
  DATA.variants.forEach(v => v.figures.forEach(f => { if(!seen.has(f.group)){ seen.add(f.group); groups.push(f.group); }}));
  const ff = document.getElementById("figfilters");

  const grow = document.createElement("div"); grow.className = "axis";
  const gnm = document.createElement("span"); gnm.className = "name"; gnm.textContent = "fields"; grow.appendChild(gnm);
  groups.forEach(g => {
    const b = document.createElement("button"); b.className = "v on"; b.textContent = g;
    b.onclick = () => { b.classList.toggle("on"); if(figOff.has(g)) figOff.delete(g); else figOff.add(g); render(); };
    grow.appendChild(b);
  });
  ff.appendChild(grow);

  const krow = document.createElement("div"); krow.className = "axis";
  const knm = document.createElement("span"); knm.className = "name"; knm.textContent = "show"; krow.appendChild(knm);
  [["static", "static frame"], ["gif", "gif"]].forEach(pair => {
    const k = pair[0], label = pair[1];
    const b = document.createElement("button"); b.className = "v on"; b.textContent = label;
    b.onclick = () => { b.classList.toggle("on"); if(kindOff.has(k)) kindOff.delete(k); else kindOff.add(k); render(); };
    krow.appendChild(b);
  });
  ff.appendChild(krow);
}

function matches(v){
  return DATA.order.every(k => sel[k].size === 0 || sel[k].has(v.params[k]));
}

// A value on axis k is "available" if some variant has it while satisfying the
// selections on every OTHER axis. Unavailable values are greyed out and disabled,
// except ones already selected (kept clickable so they can be toggled off).
function updateAvail(){
  DATA.order.forEach(k => axisBtns[k].forEach(({val, btn}) => {
    const ok = DATA.variants.some(v => v.params[k] === val &&
      DATA.order.every(j => j === k || sel[j].size === 0 || sel[j].has(v.params[j])));
    btn.disabled = !ok && !sel[k].has(val);
  }));
}

function tile(v, f){
  const cap = v.name + " — " + f.label;
  if(f.kind === "gif"){
    const gif = enc(v.name+"/"+f.file);
    const inner = f.poster
      ? '<img loading="lazy" src="'+enc(v.name+"/"+f.poster)+'" data-gif="'+gif+'" data-cap="'+cap+'" alt="'+f.label+'">'
      : '<div class="noposter" data-gif="'+gif+'" data-cap="'+cap+'">&#9654; load gif</div>';
    return '<figure><span class="badge">&#9654; gif</span>'+inner+'<figcaption>'+f.label+'</figcaption></figure>';
  }
  return '<figure><img loading="lazy" src="'+enc(v.name+"/"+f.file)+'" data-cap="'+cap+'" alt="'+f.label+'">'+
         '<figcaption>'+f.label+'</figcaption></figure>';
}

function card(v){
  const figs = v.figures.filter(f => !figOff.has(f.group) && !kindOff.has(f.kind === "gif" ? "gif" : "static"));
  const pstr = DATA.order.map(k => '<span class="k">'+k+'</span> '+(v.params[k]??'-')).join(" &nbsp; ");
  const tiles = figs.map(f => tile(v, f)).join("");
  return '<div class="card"><h2>'+pstr+'</h2><div class="grid">'+
         (tiles||'<span class="empty">no figures for selected types</span>')+'</div></div>';
}

function render(){
  const m = DATA.variants.filter(matches);
  document.getElementById("count").textContent = "showing " + m.length + " / " + DATA.variants.length + " variants";
  document.getElementById("results").innerHTML =
    m.length ? m.map(card).join("") : '<div class="empty">no variants match the selected parameters</div>';
  updateAvail();
}

document.getElementById("reset").onclick = () => {
  Object.values(sel).forEach(s => s.clear());
  figOff.clear(); kindOff.clear();
  document.querySelectorAll("#filters button.v").forEach(b => b.classList.remove("on"));
  document.querySelectorAll("#figfilters button.v").forEach(b => b.classList.add("on"));
  render();
};

const lb = document.getElementById("lb"), lbimg = lb.querySelector("img"), lbcap = lb.querySelector(".cap");
document.getElementById("results").addEventListener("click", e => {
  const t = e.target;
  const src = t.dataset.gif || (t.tagName === "IMG" ? t.getAttribute("src") : null);
  if(!src) return;
  const cap = t.dataset.cap || "";
  lbcap.textContent = t.dataset.gif ? cap + "  (loading gif…)" : cap;
  lbimg.onload = () => { lbcap.textContent = cap; };
  lbimg.src = src;            // gif fetched here, only on click
  lb.style.display = "flex";
});
lb.onclick = () => { lb.style.display = "none"; lbimg.src = ""; };
document.addEventListener("keydown", e => { if(e.key === "Escape"){ lb.style.display = "none"; lbimg.src = ""; }});

buildFilters();
render();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
