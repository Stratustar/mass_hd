#!/usr/bin/env python
"""Acceptance checks for the confluent-wet model (step 1 + chi/m dynamics).

Every check is a property the model MUST have if the implementation is right, stated
so that a number can be compared against a tolerance rather than a picture eyeballed.

  passive     zeta = 0: u -> 0, n -> rho, both pressures -> 0, q^2 -> 1, chi frozen
  advect      D_t chi = 0: chi is a pure Lagrangian tracer -- range, mean, spectrum
  turbulence  active turbulence + the trace/traceless identity p_LB - P == .5 CC (1-q^2)^2
  defectpair  the +1/2 defect self-propels
  loop        the P -> m -> chi cascade actually moves

Run:  conda run --no-capture-output -n env1 python plot/python/confluent_wet_validate.py <case_root>
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from archive.archive import loadarchive


def g(fr, name):
    lx, ly = fr.parameters["LX"], fr.parameters["LY"]
    return np.array(getattr(fr, name)).reshape((lx, ly))


def q2(fr):
    return g(fr, "QQxx")**2 + g(fr, "QQyx")**2


def vrms(fr, xf="ux_mat", yf="uy_mat"):
    return float(np.sqrt(np.mean(g(fr, xf)**2 + g(fr, yf)**2)))


def line(ok, label, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")
    return ok


def report(name, ar):
    p = ar.parameters
    nf = ar._nframes
    print(f"\n=== {name}  ({nf} frames, {p['LX']}x{p['LY']}, zeta={p['zeta']}, "
          f"CC={p['CC']}, LL={p['LL']}, tau_chi={p['tau_chi']}, tau_m={p['tau_m']})")
    xi_N = np.sqrt(p["LL"]/(2*p["CC"]))
    la = np.sqrt(p["LL"]/p["zeta"]) if p["zeta"] > 0 else np.inf
    print(f"  xi_N = {xi_N:.2f}   l_a = {la:.1f}   eta = {p['rho']/3*(p['tau']-.5):.2f}")
    return p, nf


def check_trace_split(ar, nf):
    """p_LB - P must equal the stored sigma_bulk at every node.

    This is the direct test that the code's stress really is split into an exact trace
    (sigma_bulk) plus an exactly traceless remainder, which is what makes P = -1/2 Tr(Pi)
    meaningful at all. Compare against the STORED sigma_bulk, not one recomputed from the
    stored Q: a frame's Q, chi and m are one predictor-corrector sub-step newer than its
    stress/pressure/velocity fields (Step() ends with UpdateFields, after the last
    UpdateQuantities), so recomputing gives an O(dQ) offset that has nothing to do with the
    trace split. The residual floor is the ~6 significant digits of the JSON writer.
    """
    fr = ar[nf-1]
    lhs = g(fr, "pressure_lb") - g(fr, "pressure")
    rhs = g(fr, "sigma_bulk")
    err = float(np.max(np.abs(lhs - rhs)))
    scale = max(float(np.max(np.abs(g(fr, "pressure_lb")))), 1e-30)
    tol = max(1e-12, 1e-5*scale)          # 6 significant digits of the writer
    return line(err < tol, "trace/traceless split (p_LB - P == sigma_bulk)",
                f"max err = {err:.3e}  (writer floor ~{tol:.1e})")


def check_material_velocity(ar, nf):
    """u_mat must differ from the bare LB moment, and by the amount the theory says.

    This is not a correctness bound, it is a magnitude report: if the two are identical
    the correction was never applied."""
    fr = ar[nf-1]
    du = np.hypot(g(fr, "ux_mat") - g(fr, "ux"), g(fr, "uy_mat") - g(fr, "uy"))
    u = np.hypot(g(fr, "ux_mat"), g(fr, "uy_mat"))
    urms = float(np.sqrt(np.mean(u**2)))
    print(f"  [info] |u_mat - u_code|: mean {du.mean():.3e}  max {du.max():.3e}"
          f"   (u_rms = {urms:.3e}, max ratio {du.max()/max(urms,1e-30):.2f})")


def case_passive(ar, nf):
    ok = True
    fr = ar[nf-1]
    p = ar.parameters
    ok &= line(vrms(fr) < 1e-6, "u -> 0", f"v_rms = {vrms(fr):.3e}")
    dn = np.abs(g(fr, "n") - p["rho"]).max()
    ok &= line(dn < 1e-6, "n -> rho", f"max|n-rho| = {dn:.3e}")
    for f in ("pressure", "pressure_lb"):
        v = np.abs(g(fr, f)).max()
        ok &= line(v < 1e-6, f"{f} -> 0", f"max = {v:.3e}")
    dq = np.abs(q2(fr) - 1).max()
    ok &= line(dq < 1e-3, "q^2 -> 1 (all defects annihilated)", f"max|q^2-1| = {dq:.3e}")
    dchi = np.abs(g(fr, "chi") - p["chi0"]).max()
    ok &= line(dchi < 1e-14, "chi frozen at chi0 (no spurious transport)",
               f"max|chi-chi0| = {dchi:.3e}")
    return ok


def case_advect(ar, nf):
    """D_t chi = 0 exactly, so chi is a pure Lagrangian tracer."""
    ok = True
    f0, fN = ar[0], ar[nf-1]
    c0, cN = g(f0, "chi"), g(fN, "chi")

    # 1. the clamp must never have had to fire
    eps = min(cN.min(), 1.0 - cN.max())
    ok &= line(cN.min() > -1e-12 and cN.max() < 1 + 1e-12, "chi stays in [0,1]",
               f"range [{cN.min():.6f}, {cN.max():.6f}], margin to the clamp {eps:.3e}")
    at_clamp = int(np.sum((cN <= 1e-12) | (cN >= 1 - 1e-12)))
    ok &= line(at_clamp == 0, "clamp never engaged",
               f"{at_clamp} nodes sitting on 0 or 1")

    # 2. the mean of a tracer in incompressible flow is conserved
    drift = abs(cN.mean() - c0.mean())
    ok &= line(drift < 1e-6, "<chi> conserved",
               f"|<chi>_N - <chi>_0| = {drift:.3e}  ({c0.mean():.8f} -> {cN.mean():.8f})")

    # 3. variance may only fall through the physical cascade; report the trajectory
    var = [float(np.var(g(ar[i], "chi"))) for i in range(0, nf, max(1, nf//8))]
    print("  [info] var(chi) over time: " + " ".join(f"{v:.4f}" for v in var))

    # 4. high-k content: if centred advection is dispersing, power piles up at k_max
    F = np.abs(np.fft.rfft2(cN - cN.mean()))**2
    kx = np.fft.fftfreq(cN.shape[0])[:, None]
    ky = np.fft.rfftfreq(cN.shape[1])[None, :]
    k = np.hypot(kx, ky)
    hi = float(F[k > 0.35].sum()/max(F.sum(), 1e-30))
    ok &= line(hi < 0.05, "no grid-scale pile-up in chi",
               f"fraction of chi power at k > 0.35 (wavelength < 3 cells) = {hi:.4f}")
    return ok


def case_turbulence(ar, nf):
    ok = True
    p = ar.parameters
    ok &= check_trace_split(ar, nf)
    check_material_velocity(ar, nf)
    v = [vrms(ar[i]) for i in range(nf)]
    late = np.array(v[nf//2:])
    ok &= line(late.mean() > 1e-4 and v[-1] > 0.3*max(v),
               "flow is active (v_rms saturates, not decaying)",
               f"v_rms: t0 {v[0]:.3e} -> mid {v[nf//2]:.3e} -> end {v[-1]:.3e} "
               f"(peak {max(v):.3e})")
    ma = max(v)/np.sqrt(1/3.)
    ok &= line(ma < 0.05, "Mach number small enough for the weakly-compressible limit",
               f"Ma_peak = u_rms/cs = {ma:.4f}")
    melt = [float(np.mean(q2(ar[i]) < 0.5)) for i in range(0, nf, max(1, nf//20))]
    ok &= line(max(melt) > 0, "defects appear at some point",
               f"peak melted fraction (q^2 < 0.5) over the run = {max(melt):.4f}, "
               f"final = {melt[-1]:.4f}")
    fr = ar[nf-1]
    P = g(fr, "pressure")
    print(f"  [info] P: mean {P.mean():+.3e}  std {P.std():.3e}  "
          f"[{P.min():+.3e}, {P.max():+.3e}]   (zeta = {p['zeta']})")
    print(f"  [info] suggested pmem = median(P) = {np.median(P):+.4e}, "
          f"pmem-width ~ IQR/2 = {(np.percentile(P,75)-np.percentile(P,25))/2:.4e}")
    return ok


def case_defectpair(ar, nf):
    ok = True
    check_material_velocity(ar, nf)
    v = [vrms(ar[i]) for i in range(nf)]
    ok &= line(max(v) > 1e-5, "the pair drives a flow",
               f"v_rms max = {max(v):.3e} (t0 {v[0]:.3e}, end {v[-1]:.3e})")
    qmin = [float(q2(ar[i]).min()) for i in range(0, nf, max(1, nf//10))]
    ok &= line(min(qmin) < 0.1, "the defect cores melt (xi_N resolves them)",
               f"min q^2 over the run = {min(qmin):.5f}")
    return ok


def case_loop(ar, nf):
    ok = True
    p = ar.parameters
    chi = [float(np.mean(g(ar[i], "chi"))) for i in range(nf)]
    mm = [float(np.mean(g(ar[i], "m"))) for i in range(nf)]
    PP = [float(np.mean(g(ar[i], "pressure"))) for i in range(nf)]
    ok &= line(abs(mm[-1] - p["m0"]) > 1e-6, "m responds to P",
               f"<m>: {mm[0]:.4f} -> {mm[-1]:.4f}")
    ok &= line(abs(chi[-1] - p["chi0"]) > 1e-6, "chi responds to m",
               f"<chi>: {chi[0]:.4f} -> {chi[-1]:.4f}")
    cN, mN = g(ar[nf-1], "chi"), g(ar[nf-1], "m")
    ok &= line(cN.min() >= 0 and cN.max() <= 1 and mN.min() >= 0 and mN.max() <= 1,
               "chi, m stay in [0,1]",
               f"chi [{cN.min():.4f},{cN.max():.4f}]  m [{mN.min():.4f},{mN.max():.4f}]")
    step = max(1, nf//8)
    print("  [info] <P>, <m>, <chi> over time:")
    for i in range(0, nf, step):
        print(f"         t={i*p['ninfo']:6d}  P={PP[i]:+.4e}  m={mm[i]:.4f}  chi={chi[i]:.4f}")
    return ok


CASES = {"passive": case_passive, "advect": case_advect, "turbulence": case_turbulence,
         "defectpair": case_defectpair, "loop": case_loop}

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else \
        "/Users/helu/Research/ActiveAdaptiveMatter/Code/localtest/cases/20260826/confluent_wet"
    allok = True
    for name, fn in CASES.items():
        path = os.path.join(root, name)
        if not os.path.isdir(path):
            print(f"\n=== {name}: MISSING ({path})")
            allok = False
            continue
        ar = loadarchive(path)
        _, nf = report(name, ar)
        try:
            allok &= bool(fn(ar, nf))
        except Exception as e:
            print(f"  [FAIL] exception: {type(e).__name__}: {e}")
            allok = False
    print("\n" + ("ALL CHECKS PASSED" if allok else "SOME CHECKS FAILED"))
    sys.exit(0 if allok else 1)
