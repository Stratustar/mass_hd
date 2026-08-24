"""Phase-separation diagnostics for the phenotype field of a confluent-memory scan.

A correlation length says how BIG the phenotype structures are; it does not say whether
they are two coexisting phases. This computes the three quantities that do:

  Sigma = var(chi) / (<chi>(1 - <chi>))   1 if chi is binary, 0 if uniform
  f_lo, f_hi                              area fractions with chi < 0.25 and chi > 0.75
  l_chichi(t)                             still growing = true coarsening; flat = pattern

plus the chi histogram over the steady window, since bimodality is the decisive signature.

Usage:  cm_phase.py <study> [<study> ...] --out results.json [--steady 30000] [--L 400]
where <study> is a path under $MASS_SCRATCH/cases, e.g. 20260822/cm_regime.
"""
import argparse, json, os
import numpy as np
from cm_common import (SCRATCH, frame_count, iter_frames, Radial, parse_variant, ph)

BINS = np.linspace(0, 1, 51)


def analyse_case(root, radial, steady):
    if frame_count(root) < 20:
        return None
    hist = np.zeros(len(BINS) - 1)
    ts = {k: [] for k in ("t", "l", "S", "flo", "fhi", "chibar")}
    nsteady = 0
    for _, t, fr in iter_frames(root):
        chi = ph.grid(fr, "chi")
        mu, var = float(chi.mean()), float(chi.var())
        c = radial.mean(radial.corr([chi]))
        c = c / c[0] if c[0] > 0 else c * np.nan
        ts["t"].append(t)
        ts["l"].append(radial.length(c))
        ts["S"].append(var / max(mu * (1 - mu), 1e-12))
        ts["flo"].append(float((chi < 0.25).mean()))
        ts["fhi"].append(float((chi > 0.75).mean()))
        ts["chibar"].append(mu)
        if t >= steady:
            hist += np.histogram(chi.ravel(), bins=BINS)[0]
            nsteady += 1
    if nsteady == 0:
        return None
    hist /= hist.sum()
    st = [j for j, t in enumerate(ts["t"]) if t >= steady]
    out = {"hist": hist.tolist(), "bins": BINS.tolist(), "nsteady": nsteady}
    for k, v in ts.items():
        out["ts_" + k] = v
        if k != "t":
            out["mean_" + k] = float(np.nanmean([v[j] for j in st]))
    # coarsening: fit l ~ t^n over the second half, and compare first and last steady values
    tt = np.array([ts["t"][j] for j in st], float)
    ll = np.array([ts["l"][j] for j in st], float)
    ok = np.isfinite(ll) & (tt > 0) & (ll > 0)
    out["coarsen_exponent"] = (float(np.polyfit(np.log(tt[ok]), np.log(ll[ok]), 1)[0])
                               if ok.sum() > 5 else float("nan"))
    out["l_first_last"] = [float(ll[ok][0]), float(ll[ok][-1])] if ok.sum() else [np.nan, np.nan]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--steady", type=float, default=None,
                    help="steady-window start; default = half the run of each case")
    ap.add_argument("--L", type=int, default=400)
    args = ap.parse_args()

    radial = Radial(args.L)
    results = []
    for study in args.studies:
        base = os.path.join(SCRATCH, "cases", study)
        for case in sorted(os.listdir(base)):
            root = os.path.join(base, case)
            if not os.path.isdir(root):
                continue
            steady = args.steady
            if steady is None:
                oa_n = frame_count(root)
                _, tmax, _ = list(iter_frames(root))[-1] if oa_n else (0, 0, None)
                steady = tmax / 2.
            r = analyse_case(root, radial, steady)
            if r is None:
                print(f"skip {study}/{case} (too few frames)", flush=True)
                continue
            r.update(study=study, case=case, steady=steady, **parse_variant(case))
            results.append(r)
            print(f"done {study}/{case}  Sigma={r['mean_S']:.3f} "
                  f"l={r['mean_l']:.1f} f_lo={r['mean_flo']:.3f} f_hi={r['mean_fhi']:.3f}",
                  flush=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(results, open(args.out, "w"))
    print(f"\nwrote {len(results)} cases -> {args.out}")


if __name__ == "__main__":
    main()
