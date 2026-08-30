#!/usr/bin/env python3
"""Run the confluent-wet scan analysis for ONE case, and write it as a part file.

cw_scan.py loops over a whole study in a single process, which is right at L=400 -- a few
minutes a case -- and hopeless at L=800, where reading the 50-frame analysis window costs
~19 min per case and a 248-case scan would take 79 hours in one job. Splitting it one case
per array task turns that into one 19-minute task run 20-wide.

Takes the array runner's contract (`<case_input_dir> <case_output_dir>`) so it can be chained
straight behind the simulation array, and writes `<case_output_dir>/part.json`. cw_merge_scan.py
then concatenates the parts into the single cw_scan.json that cw_figs and cw_html expect.

Usage: cw_scan_one.py <case_input_dir> <case_output_dir> [--nlast 50] [--stride 3]
                      [--maxlag-frac 0.4]
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cw_scan import analyse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputdir")
    ap.add_argument("outdir")
    ap.add_argument("--nlast", type=int, default=50)
    ap.add_argument("--stride", type=int, default=3)
    ap.add_argument("--maxlag-frac", type=float, default=0.4)
    a = ap.parse_args()

    root = a.inputdir.rstrip("/")
    case = os.path.basename(root)
    r = analyse(root, a.nlast, a.stride, a.maxlag_frac)
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "part.json")
    with open(out, "w") as f:
        json.dump(r, f)
    print(f"{case}: A_eff={r['A_eff']:.2f} <chi>={r['chi_bar']:.3f} "
          f"t_eddy={r['t_eddy']:.0f} L_chi={r['L_chi']:.1f} "
          f"lag(P->m)={r['lag_Pm_k0']:.0f}/{r['tau_m_set']:.0f} "
          f"lag(m->chi)={r['lag_mchi_k0']:.0f}/{r['tau_chi_set']:.0f}", flush=True)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
