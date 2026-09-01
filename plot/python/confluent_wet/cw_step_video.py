#!/usr/bin/env python3
"""Build the video dashboard for the step-switch tau_m campaign.

    cw_step_video.py <b1_results> <b2_results> --videos DIR --calib calib.json --out DIR

WHY A LOCAL PAGE AND NOT AN ARTIFACT. The campaign wrote 400 mp4s totalling 25 GB, and an
Artifact cannot host video: its CSP blocks every external host, and inlining as data URIs
runs into a 16 MB page cap several hundred times over. So this is a plain HTML file that
sits beside the (transcoded) clips and is opened from disk.

WHAT IT SHOWS. The chi field for all 100 production runs, laid out exactly like the fate
map -- rows are tau_m, columns are the four initial conditions, one section per tau_chi
rule -- because the whole point of these videos is COMPARISON, and the comparison that
matters is across a row: four starts at one memory time. Each row therefore gets its own
sync-play control. A handful of runs additionally carry m, P and |u|, for the cases where
the question is what the flow is doing rather than what chi settled to.

Videos are lazily attached: 100 <video> elements with real sources would have the browser
fetching a quarter of a gigabyte on load. An IntersectionObserver sets src when a card
comes near the viewport, and preload stays at metadata until something asks to play.
"""
import argparse
import glob
import html
import json
import os
import re

FIELDS = [("chi", "χ 表型", "0 → 1"), ("m", "m 记忆", "0 → 1"),
          ("P", "P 压强", "±3 σ_P"), ("u", "|u| 速度", "0 → 3 u_rms")]
STARTS = [("a", "χ ≡ 0", "活性初值"), ("b", "χ ≡ 1", "被动初值"),
          ("c1", "二值噪声 1", "chi-seed 1"), ("c2", "二值噪声 2", "chi-seed 2")]
FATE_CLASS = {"active": "f-act", "passive": "f-pas",
              "undecided": "f-und", "mixed": "f-mix"}
FATE_LABEL = {"active": "活性", "passive": "被动", "undecided": "未定", "mixed": "混合"}


def load(d):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        case = os.path.basename(os.path.dirname(p))
        m = re.match(r"tm(.+)_([abc]\d?)$", case)
        if not m:
            continue
        with open(p) as fh:
            out[(float(m.group(1).replace("p", ".")), m.group(2))] = (case, json.load(fh))
    return out


def card(grp, case, start_key, d, vdir, field="chi"):
    """One video card: the clip plus the numbers that say what it settled to."""
    f = d.get("fate") or {}
    s = d.get("settled") or {}
    fl = d["flow"]
    u = d["params"]
    label = f.get("label", "?")
    chi = f.get("chi_tail", s.get("chi_mean_tail", float("nan")))
    dep = f.get("departure_steps", float("nan"))
    life = (f"{dep/u['tau_m']:.1f} τ_m" if dep == dep else "从未离开")
    src = f"clips/{grp}_{case}_{field}.mp4"
    exists = os.path.exists(os.path.join(vdir, f"{grp}_{case}_{field}.mp4"))
    if not exists:
        return (f'<figure class="card missing"><div class="ph">缺 {field}</div>'
                f'<figcaption><div class="ct">{html.escape(case)}</div></figcaption></figure>')
    return f"""<figure class="card {FATE_CLASS.get(label,'f-mix')}">
  <video data-src="{src}" muted playsinline loop preload="none"></video>
  <figcaption>
    <div class="ct"><span>{html.escape(start_key)}</span><span class="badge">{FATE_LABEL.get(label,label)}</span></div>
    <div class="stats"><span class="st"><b>⟨χ⟩</b> {chi:.3f}</span>
      <span class="st"><b>寿命</b> {life}</span>
      <span class="st"><b>N_def</b> {fl['N_def']:.0f}</span>
      <span class="st"><b>melt</b> {fl['melted_frac']:.0%}</span></div>
  </figcaption>
</figure>"""


def section(grp, title, subtitle, data, vdir):
    gs = sorted(set(g for g, _ in data))
    rows = []
    for g in gs:
        cells = []
        for key, lab, _sub in STARTS:
            if (g, key) not in data:
                cells.append('<figure class="card missing"><div class="ph">—</div>'
                             '<figcaption><div class="ct">未跑</div></figcaption></figure>')
                continue
            case, d = data[(g, key)]
            cells.append(card(grp, case, lab, d, vdir))
        _case, d0 = data[(gs[0], "a")]
        tm_steps = data[(g, "a")][1]["params"]["tau_m"]
        rows.append(f"""<div class="row" data-row>
  <div class="rlab"><div class="rt">τ_m = {g:g} τ_c</div>
    <div class="rs">{tm_steps:.0f} 步</div>
    <button class="sync" type="button">▶ 同步播放这一行</button></div>
  <div class="grid4">{''.join(cells)}</div>
</div>""")
    return f"""<section>
  <div class="shead"><h2>{title}</h2><span class="cs">{subtitle}</span></div>
  {''.join(rows)}
</section>"""


def deep_section(picks, b1, b2, vdir):
    blocks = []
    for grp, case in picks:
        data = b1 if grp == "b1" else b2
        hit = [(k, v) for k, v in data.items() if v[0] == case]
        if not hit:
            continue
        (g, key), (_c, d) = hit[0]
        cells = []
        for fkey, flab, fscale in FIELDS:
            if not os.path.exists(os.path.join(vdir, f"{grp}_{case}_{fkey}.mp4")):
                continue
            cells.append(f"""<figure class="card">
  <video data-src="clips/{grp}_{case}_{fkey}.mp4" muted playsinline loop preload="none"></video>
  <figcaption><div class="ct">{flab}</div><div class="cs">{fscale}</div></figcaption>
</figure>""")
        if not cells:
            continue
        f = d.get("fate") or {}
        start = dict((k, l) for k, l, _ in STARTS).get(key, key)
        blocks.append(f"""<div class="row" data-row>
  <div class="rlab"><div class="rt">{grp.upper()} · τ_m = {g:g} τ_c</div>
    <div class="rs">{start} → {FATE_LABEL.get(f.get('label','?'),'?')} ({f.get('chi_tail',float('nan')):.3f})</div>
    <button class="sync" type="button">▶ 同步播放这一行</button></div>
  <div class="grid4">{''.join(cells)}</div>
</div>""")
    return f"""<section>
  <div class="shead"><h2>四场对照</h2><span class="cs">同一跑的 χ / m / P / |u|，选出来的十个跑</span></div>
  {''.join(blocks)}
</section>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("b1")
    ap.add_argument("b2")
    ap.add_argument("--videos", required=True, help="dir holding the transcoded clips")
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.calib) as fh:
        cal = json.load(fh)
    b1, b2 = load(a.b1), load(a.b2)
    os.makedirs(a.out, exist_ok=True)

    picks = [("b1", "tm1_a"), ("b1", "tm4p7_a"), ("b1", "tm10_a"), ("b1", "tm22_a"),
             ("b1", "tm30_a"), ("b1", "tm30_b"), ("b1", "tm30_c1"),
             ("b2", "tm22_a"), ("b2", "tm22_b"), ("b2", "tm30_a")]

    tau_x = cal.get("tau_x_coexistence") or cal.get("tau_x")
    body = f"""<header>
  <div class="eyebrow">confluent-wet · 2026-09-01 · 阶跃版力学记忆 · 帧距 0.4 τ_c</div>
  <h1>τ_m 扫描 · 视频对照台</h1>
  <p class="lede">每一行是一个 τ_m 下的四种初值，<b>可以整行同步播放</b>——这些视频的用处就是横向对照。
  点画面可放大。色标全场统一：χ 与 m 在 [0,1]，P 在 ±3σ_P，|u| 在 [0, 3u_rms]。
  <b>χ 是 coolwarm：蓝 = 0 = 活性，红 = 1 = 被动</b>，命运图用的是同一套方向。</p>
  <div class="scales">
    <div class="scale"><div class="k">τ_c(ζ)</div><div class="v">{cal['tau_c']:.0f} 步</div></div>
    <div class="scale"><div class="k">mc</div><div class="v">{cal['mc']:.4f}</div></div>
    <div class="scale"><div class="k">窗口 (f_floor, f_top)</div>
      <div class="v">({cal['f_floor']:.4f}, {cal['f_top']:.4f})</div></div>
    <div class="scale"><div class="k">平均场阈值</div><div class="v">{tau_x:.1f} τ_c</div></div>
  </div>
  <div class="scales">
    <div class="scale"><div class="k">χ · coolwarm</div>
      <div class="v">0 活性 → 1 被动</div>
      <div class="ramp" style="background:linear-gradient(90deg,#3b4cc0,#dddddd,#b40426)"></div></div>
    <div class="scale"><div class="k">m · viridis</div><div class="v">0 → 1</div>
      <div class="ramp" style="background:linear-gradient(90deg,#440154,#21918c,#fde725)"></div></div>
    <div class="scale"><div class="k">P · RdBu_r</div><div class="v">±{3*cal['sigma_P']:.4f}</div>
      <div class="ramp" style="background:linear-gradient(90deg,#053061,#f7f7f7,#67001f)"></div></div>
    <div class="scale"><div class="k">|u| · magma</div><div class="v">0 → {3*cal['u_rms']:.4f}</div>
      <div class="ramp" style="background:linear-gradient(90deg,#000004,#b5367a,#fcfdbf)"></div></div>
  </div>
  <div class="legend">
    <span class="lg f-act">活性</span><span class="lg f-pas">被动</span>
    <span class="lg f-und">未定（仍在漂移）</span><span class="lg f-mix">混合</span>
  </div>
</header>
{section('b1', 'B1 · τ_χ = 0.3 τ_c', '13 个 τ_m × 4 种初值，χ 场', b1, a.videos)}
{section('b2', 'B2 · τ_χ = τ_m', '12 个 τ_m × 4 种初值，χ 场', b2, a.videos)}
{deep_section(picks, b1, b2, a.videos)}"""

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "cw_step_video.css")) as fh:
        css = fh.read()
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "cw_step_video.js")) as fh:
        js = fh.read()
    page = (f"<!doctype html>\n<html lang=\"zh\"><head><meta charset=\"utf-8\">"
            f"<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<title>τ_m 扫描 · 视频对照台</title>\n<style>{css}</style></head>\n"
            f"<body><div class=\"wrap\">{body}</div>"
            f"<dialog id=\"zoom\"><video controls loop autoplay muted playsinline></video>"
            f"<div class=\"zt\"></div></dialog>\n<script>{js}</script></body></html>\n")
    out = os.path.join(a.out, "index.html")
    with open(out, "w") as fh:
        fh.write(page)
    n = len(glob.glob(os.path.join(a.videos, "*.mp4")))
    print(f"wrote {out} ({len(page)/1024:.0f} KB) against {n} clips")


if __name__ == "__main__":
    main()
