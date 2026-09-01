#!/usr/bin/env python3
"""Video dashboard for the symmetry test: three operating points side by side.

    cw_step_sym_video.py <b1_results> <sym_results> --videos DIR --calib calib.json --out DIR

WHY A SEPARATE PAGE from cw_step_video.py. That one lays the 100 production runs out like
the fate map -- rows are tau_m, columns are the four initial conditions -- because the
question there was "what does each start do". Here the question is different and so is the
comparison that answers it: at ONE tau_m, what do the two pure starts do at the campaign's
mc against the two symmetric operating points. So a row is a tau_m and the six columns are
(base, s1, s2) x (chi = 0, chi = 1), and reading across it shows the whole result --
the campaign's active arm collapsing where both symmetric points hold it.

Clips are shared with the main dashboard's directory; the base columns are the very files
that page already uses, not copies.
"""
import argparse
import glob
import json
import os
import re

VARIANTS = [
    ("base", "b1", "原 mc = 0.2610", "cusp 选出的阈值"),
    ("s1", "sym_s1", "s1 · mc = 0.2100", "pmem 不变，阈值居中"),
    ("s2", "sym_s2", "s2 · mc = 0.1585", "pmem = 0.6 σ_P"),
]
ARMS = [("a", "χ ≡ 0"), ("b", "χ ≡ 1")]
TAU_M = [10, 15, 22, 30]
FIELDS = [("chi", "χ 表型", "coolwarm 0→1"), ("m", "m 记忆", "viridis 0→1"),
          ("P", "P 压强", "±3 σ_P"), ("u", "|u| 速度", "0 → 3 u_rms")]
FATE_CLASS = {"active": "f-act", "passive": "f-pas",
              "undecided": "f-und", "mixed": "f-mix"}
FATE_LABEL = {"active": "活性", "passive": "被动", "undecided": "未定", "mixed": "混合"}


def load_base(d):
    """{(tau_m, arm): part} from the campaign's B1 tree."""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        c = os.path.basename(os.path.dirname(p))
        m = re.match(r"tm(.+)_([ab])$", c)
        if not m:
            continue
        with open(p) as fh:
            out[(float(m.group(1).replace("p", ".")), m.group(2))] = (c, json.load(fh))
    return out


def load_sym(d):
    """{(variant, tau_m, arm): part} from the symmetry tree."""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        c = os.path.basename(os.path.dirname(p))
        m = re.match(r"(s\d)_tm(.+)_([ab])$", c)
        if not m:
            continue
        with open(p) as fh:
            out[(m.group(1), float(m.group(2).replace("p", ".")), m.group(3))] = \
                (c, json.load(fh))
    return out


def margin(d, arm):
    """The phase's distance from the threshold, in its OWN sigma_m -- the whole point."""
    s = d.get("settled") or {}
    u = d["params"]
    m, sd, mc = s.get("m_mean_tail", float("nan")), d["flow"]["std_m"], u["mc"]
    if sd <= 0 or m != m:
        return float("nan")
    return (m - mc) / sd if arm == "a" else (mc - m) / sd


def card(clip, d, arm, vdir, title, field="chi"):
    f = d.get("fate") or {}
    s = d.get("settled") or {}
    u = d["params"]
    label = f.get("label", "?")
    chi = f.get("chi_tail", s.get("chi_mean_tail", float("nan")))
    dep = f.get("departure_steps", float("nan"))
    life = f"{dep/u['tau_m']:.1f} τ_m" if dep == dep else "从未离开"
    z = margin(d, arm)
    name = f"{clip}_{field}.mp4"
    if not os.path.exists(os.path.join(vdir, name)):
        return ('<figure class="card missing"><div class="ph">缺片段</div>'
                f'<figcaption><div class="ct"><span>{title}</span></div></figcaption></figure>')
    zs = f"{z:.2f} σ" if z == z else "—"
    return f"""<figure class="card {FATE_CLASS.get(label,'f-mix')}">
  <video data-src="clips/{name}" muted playsinline loop preload="none"></video>
  <figcaption>
    <div class="ct"><span>{title}</span><span class="badge">{FATE_LABEL.get(label,label)}</span></div>
    <div class="stats"><span class="st"><b>⟨χ⟩</b> {chi:.4f}</span>
      <span class="st"><b>余量</b> {zs}</span>
      <span class="st"><b>寿命</b> {life}</span></div>
  </figcaption>
</figure>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("b1")
    ap.add_argument("sym")
    ap.add_argument("--videos", required=True)
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.calib) as fh:
        cal = json.load(fh)
    base, sym = load_base(a.b1), load_sym(a.sym)
    os.makedirs(a.out, exist_ok=True)

    rows = []
    for g in TAU_M:
        cells = []
        for vkey, prefix, vlab, _sub in VARIANTS:
            for arm, alab in ARMS:
                if vkey == "base":
                    hit = base.get((g, arm))
                    clip = f"b1_{hit[0]}" if hit else None
                else:
                    hit = sym.get((vkey, g, arm))
                    clip = f"sym_{hit[0]}" if hit else None
                if not hit:
                    cells.append('<figure class="card missing"><div class="ph">未跑</div>'
                                 '<figcaption><div class="ct"><span>—</span></div>'
                                 '</figcaption></figure>')
                    continue
                cells.append(card(clip, hit[1], arm, a.videos,
                                  f"{vlab.split(' · ')[0]} · {alab}"))
        tm = (base.get((g, "a")) or sym.get(("s1", g, "a")))[1]["params"]["tau_m"]
        rows.append(f"""<div class="row6" data-row>
  <div class="rlab"><div class="rt">τ_m = {g:g} τ_c</div><div class="rs">{tm:.0f} 步</div>
    <button class="sync" type="button">▶ 同步播放这一行</button></div>
  <div class="grid6">{''.join(cells)}</div>
</div>""")

    deep = []
    for case, what in (("s1_tm10_a", "s1 · τ_m = 10 τ_c · χ≡0 —— 原参数在这里塌了，这里守住"),
                       ("s1_tm10_b", "s1 · τ_m = 10 τ_c · χ≡1 —— 反向塌陷：被动臂走向活性"),
                       ("s1_tm30_a", "s1 · τ_m = 30 τ_c · χ≡0 —— 对称双稳的活性相"),
                       ("s1_tm30_b", "s1 · τ_m = 30 τ_c · χ≡1 —— 对称双稳的被动相")):
        cells = []
        for fk, flab, fsc in FIELDS:
            if not os.path.exists(os.path.join(a.videos, f"sym_{case}_{fk}.mp4")):
                continue
            cells.append(f"""<figure class="card">
  <video data-src="clips/sym_{case}_{fk}.mp4" muted playsinline loop preload="none"></video>
  <figcaption><div class="ct"><span>{flab}</span></div><div class="cs">{fsc}</div></figcaption>
</figure>""")
        if cells:
            deep.append(f"""<div class="row" data-row>
  <div class="rlab"><div class="rt">{case}</div><div class="rs">{what}</div>
    <button class="sync" type="button">▶ 同步播放这一行</button></div>
  <div class="grid4">{''.join(cells)}</div>
</div>""")

    body = f"""<header>
  <div class="eyebrow">confluent-wet · 2026-09-01 · 对称化测试 · 帧距 0.4 τ_c</div>
  <h1>三个工作点，同一条 τ_m</h1>
  <p class="lede">战役参数下两个相差 2.7 倍——活性相离阈值 2.55 σ_m、被动相 6.97 σ_m，
  所以衰变的永远是活性相。把 (m − mc)/σ 在两侧配平得到 <b>mc = 0.21</b>（s1）。
  每一行是一个 τ_m，六列是三个工作点 × 两种纯初值，<b>整行同步播放</b>就能看到
  左边两格塌掉而右边四格守住。χ 是 coolwarm：<b>蓝 = 活性，红 = 被动</b>。</p>
  <div class="scales">
    <div class="scale"><div class="k">τ_c(ζ)</div><div class="v">{cal['tau_c']:.0f} 步</div></div>
    <div class="scale"><div class="k">原 mc</div><div class="v">0.2610</div>
      <div class="v" style="font-size:11px">2.55 / 6.97 σ</div></div>
    <div class="scale"><div class="k">s1 mc</div><div class="v">0.2100</div>
      <div class="v" style="font-size:11px">5.27 / 5.26 σ @ 30 τ_c</div></div>
    <div class="scale"><div class="k">s2 mc</div><div class="v">0.1585</div>
      <div class="v" style="font-size:11px">pmem 0.6 σ_P</div></div>
    <div class="scale"><div class="k">共存阈值</div><div class="v">22–30 → 10–15 τ_c</div></div>
  </div>
  <div class="legend">
    <span class="lg f-act">活性</span><span class="lg f-pas">被动</span>
    <span class="lg f-und">未定</span><span class="lg f-mix">混合</span>
  </div>
</header>
<section>
  <div class="shead"><h2>三点对照 · χ 场</h2>
    <span class="cs">列：原 / s1 / s2 × χ≡0 / χ≡1　·　「余量」是该相到自己阈值的距离，以它自己的 σ_m 计</span></div>
  {''.join(rows)}
</section>
<section>
  <div class="shead"><h2>s1 四场对照</h2><span class="cs">χ / m / P / |u|，推荐工作点的四个关键跑</span></div>
  {''.join(deep)}
</section>"""

    here = os.path.dirname(os.path.abspath(__file__))
    css = open(os.path.join(here, "cw_step_video.css")).read()
    js = open(os.path.join(here, "cw_step_video.js")).read()
    css += """
.row6{display:grid; grid-template-columns:158px 1fr; gap:16px; align-items:start;
  padding:13px 0; border-bottom:1px solid var(--rule)}
.row6:last-child{border-bottom:none}
.grid6{display:grid; grid-template-columns:repeat(6,1fr); gap:10px}
@media (max-width:1400px){.grid6{grid-template-columns:repeat(3,1fr)}}
@media (max-width:1080px){.grid6{grid-template-columns:repeat(2,1fr)}
  .row6{grid-template-columns:1fr}}
@media (max-width:560px){.grid6{grid-template-columns:1fr}}
"""
    page = ('<!doctype html>\n<html lang="zh"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>对称化测试 · 视频对照台</title>\n<style>' + css + '</style></head>\n'
            '<body><div class="wrap">' + body + '</div>'
            '<dialog id="zoom"><video controls loop autoplay muted playsinline></video>'
            '<div class="zt"></div></dialog>\n<script>' + js + '</script></body></html>\n')
    out = os.path.join(a.out, "index.html")
    with open(out, "w") as fh:
        fh.write(page)
    print(f"wrote {out} ({len(page)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
