"""Steady-state statistics, radial correlations and spectra for a confluent-memory scan.

Per case, over the steady window: time series of N_defect, v_rms, <|E|>, <chi>, std(chi),
<m>, std(m), median P and phi_min; the isotropic correlation functions C_vv, C_QQ, C_mm,
C_chichi, C_PP and the cross-correlation C_m,chi with their correlation lengths; and the
kinetic-energy spectrum. Also a detrended power spectrum of each time series, so an
oscillatory branch would show up as a peak well inside the resolvable band.

Usage:  cm_scan.py <study> --out results.json [--steady 30000] [--L 400] [--transient 20000]
"""
import argparse, json, os
import numpy as np
from cm_common import (SCRATCH, frame_count, iter_frames, Radial, parse_variant,
                       strain_rate, n_defects, ph)


def analyse_case(root, radial, steady, transient):
    if frame_count(root) < 20:
        return None
    keys = ("t", "nd", "vrms", "E", "chibar", "chistd", "mbar", "mstd", "Pmed", "phimin", "q2")
    S = {k: [] for k in keys}
    acc = {k: np.zeros(len(radial.bins)) for k in ("vv", "QQ", "mm", "cc", "PP", "mc")}
    kb = np.arange(0, radial.L // 2 + 1)
    Ek = np.zeros(len(kb))
    ky, kx = np.meshgrid(np.fft.fftfreq(radial.L) * radial.L,
                         np.fft.fftfreq(radial.L) * radial.L, indexing="ij")
    kidx = np.digitize(np.sqrt(kx**2 + ky**2).ravel(), kb) - 1
    kcnt = np.bincount(kidx, minlength=len(kb))[:len(kb)]
    n = 0
    for _, t, fr in iter_frames(root):
        ux, uy = ph.grid(fr, "ux"), ph.grid(fr, "uy")
        qxx, qyx = ph.grid(fr, "QQxx"), ph.grid(fr, "QQyx")
        chi, m = ph.grid(fr, "chi"), ph.grid(fr, "m")
        phi, P = ph.grid(fr, "phi"), -ph.grid(fr, "sigma_bulk")
        S["t"].append(t); S["nd"].append(n_defects(qxx, qyx))
        S["vrms"].append(float(np.sqrt(np.mean(ux**2 + uy**2))))
        S["E"].append(float(strain_rate(ux, uy).mean()))
        S["chibar"].append(float(chi.mean())); S["chistd"].append(float(chi.std()))
        S["mbar"].append(float(m.mean()));     S["mstd"].append(float(m.std()))
        S["Pmed"].append(float(np.median(P))); S["phimin"].append(float(phi.min()))
        S["q2"].append(float((qxx**2 + qyx**2).mean()))
        if t < steady:
            continue
        acc["vv"] += radial.mean(radial.corr([ux, uy]))
        acc["QQ"] += radial.mean(radial.corr([qxx, qyx]))
        acc["mm"] += radial.mean(radial.corr([m]))
        acc["cc"] += radial.mean(radial.corr([chi]))
        acc["PP"] += radial.mean(radial.corr([P]))
        acc["mc"] += radial.mean(radial.corr([m], [chi]))
        e = .5 * (np.abs(np.fft.fft2(ux - ux.mean()))**2
                  + np.abs(np.fft.fft2(uy - uy.mean()))**2) / (radial.L * radial.L)**2
        s = np.bincount(kidx, weights=e.ravel(), minlength=len(kb))[:len(kb)]
        Ek += np.where(kcnt > 0, s / np.maximum(kcnt, 1), 0.)
        n += 1
    if n == 0:
        return None
    for k in acc:
        acc[k] /= n
    Ek /= n

    d = {"r": radial.bins.tolist(), "k": kb.tolist(), "Ek": Ek.tolist(), "nsteady": n}
    for k in ("vv", "QQ", "mm", "cc", "PP"):
        c = acc[k] / acc[k][0] if acc[k][0] != 0 else acc[k] * np.nan
        d["C_" + k] = c.tolist(); d["len_" + k] = radial.length(c)
    den = np.sqrt(acc["mm"][0] * acc["cc"][0])
    d["C_mc"] = (acc["mc"] / den).tolist() if den > 0 else [float("nan")] * len(radial.bins)
    d["x0_mc"] = float(acc["mc"][0] / den) if den > 0 else float("nan")

    st = [j for j, t in enumerate(S["t"]) if t >= steady]
    for k in S:
        d["ts_" + k] = S[k]
        if k != "t":
            d["mean_" + k] = float(np.mean([S[k][j] for j in st]))
    tr = [j for j, t in enumerate(S["t"]) if t >= transient]
    if len(tr) > 8:
        dt = S["t"][1] - S["t"][0]
        for k in ("vrms", "nd", "chibar"):
            y = np.array([S[k][j] for j in tr], float)
            y = y - np.polyval(np.polyfit(np.arange(len(y)), y, 2), np.arange(len(y)))
            Pw = np.abs(np.fft.rfft(y * np.hanning(len(y))))**2
            f = np.fft.rfftfreq(len(y), dt)
            i = int(np.argmax(Pw[1:]) + 1)
            d["peak_period_" + k] = float(1 / f[i]) if f[i] > 0 else float("nan")
            d["peak_power_" + k] = float(Pw[i] / np.sum(Pw[1:]))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--steady", type=float, required=True)
    ap.add_argument("--transient", type=float, default=None)
    ap.add_argument("--L", type=int, default=400)
    args = ap.parse_args()
    transient = args.transient if args.transient is not None else args.steady / 2.5

    radial = Radial(args.L)
    results = []
    for study in args.studies:
        base = os.path.join(SCRATCH, "cases", study)
        for case in sorted(os.listdir(base)):
            root = os.path.join(base, case)
            if not os.path.isdir(root):
                continue
            r = analyse_case(root, radial, args.steady, transient)
            if r is None:
                print(f"skip {study}/{case}", flush=True)
                continue
            r.update(study=study, case=case, **parse_variant(case))
            results.append(r)
            print(f"done {study}/{case}", flush=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(results, open(args.out, "w"))
    print(f"\nwrote {len(results)} cases -> {args.out}")


if __name__ == "__main__":
    main()
