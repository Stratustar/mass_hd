#!/usr/bin/env python3
"""Concatenate the per-case part.json files of a scan into one cw_scan.json.

The parts come from cw_scan_one.py, one array task each; cw_figs and cw_html both want the
single list that cw_scan.py used to produce, so this restores it. Missing parts are reported
rather than skipped silently -- a scan that quietly analysed 240 of 248 cases and said nothing
is how a hole in a parameter sweep survives to the figures.

Usage: cw_merge_scan.py <results_study_dir> [--out <path>]
"""
import argparse
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resultsdir")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    root = a.resultsdir.rstrip("/")
    cases = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    res, missing = [], []
    for c in cases:
        p = os.path.join(root, c, "scan", "part.json")
        if not os.path.exists(p):
            missing.append(c)
            continue
        try:
            with open(p) as f:
                res.append(json.load(f))
        except Exception as exc:
            missing.append(f"{c} ({type(exc).__name__})")
    out = a.out or os.path.join(root, "cw_scan.json")
    with open(out, "w") as f:
        json.dump(res, f)
    print(f"merged {len(res)} parts -> {out}")
    if missing:
        print(f"MISSING or unreadable ({len(missing)}): {', '.join(missing[:20])}"
              + (" ..." if len(missing) > 20 else ""))


if __name__ == "__main__":
    main()
