#!/usr/bin/env python3
"""Check Sum(m)=Sum(phi*chi) and Sum(phi) conservation across the frames of one or
more MASS archives. Prints a per-frame table and the total relative drift, so a
transport-conservation fix can be verified (uniform box, alpha=0, switching/Dchi off:
Sum(m) should stay constant to round-off).

Usage:
  python conservation_check.py ARCHIVE_DIR [ARCHIVE_DIR ...] [--labels L1 L2 ...]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive.archive import loadarchive


def frame_indices(ar):
    n = int((ar.nsteps - ar.nstart) / ar.ninfo) + 1
    return list(range(n))


def sums(frame):
    phi = np.array(frame.phi, dtype=float)
    chi = np.array(frame.chi, dtype=float)
    m = phi * chi
    if hasattr(frame, "m"):
        m = np.array(frame.m, dtype=float)
    return phi.sum(), m.sum(), float(chi.min()), float(chi.max())


def check(archive_dir, label):
    ar = loadarchive(archive_dir)
    print(f"\n=== {label}: {archive_dir} ===")
    print(f"{'step':>8} {'Sum(phi)':>16} {'Sum(m)':>16} "
          f"{'Sum(m)/Sum(m)0':>16} {'min chi':>10} {'max chi':>10}")
    m0 = None
    phi0 = None
    rows = []
    for i in frame_indices(ar):
        step = ar.nstart + i * ar.ninfo
        try:
            fr = ar.read_frame(i)
        except Exception as exc:  # nan blow-up etc.
            print(f"{step:>8}  <unreadable: {exc}>")
            break
        sp, sm, cmin, cmax = sums(fr)
        if m0 is None:
            m0, phi0 = sm, sp
        ratio = sm / m0 if m0 else float("nan")
        rows.append((step, sp, sm, ratio, cmin, cmax))
        print(f"{step:>8} {sp:>16.9e} {sm:>16.9e} {ratio:>16.12f} "
              f"{cmin:>10.4f} {cmax:>10.4f}")
    if rows:
        sm_drift = max(abs(r[3] - 1.0) for r in rows)
        sp_drift = max(abs(r[1] / phi0 - 1.0) for r in rows) if phi0 else float("nan")
        cmax_over = max(r[5] for r in rows) - 1.0
        cmin_under = -min(r[4] for r in rows)
        print(f"  --> max |Sum(m)/Sum(m)0 - 1| = {sm_drift:.3e}")
        print(f"  --> max |Sum(phi)/Sum(phi)0 - 1| = {sp_drift:.3e}")
        print(f"  --> max(chi)-1 = {cmax_over:.3e}   -min(chi) = {cmin_under:.3e}")
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("archives", nargs="+")
    p.add_argument("--labels", nargs="*", default=None)
    args = p.parse_args()
    labels = args.labels or [os.path.basename(os.path.normpath(a)) for a in args.archives]
    if len(labels) != len(args.archives):
        p.error("number of --labels must match number of archives")
    for a, l in zip(args.archives, labels):
        check(a, l)


if __name__ == "__main__":
    main()
