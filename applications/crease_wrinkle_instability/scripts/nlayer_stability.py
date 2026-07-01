#!/usr/bin/env python3
r"""
General N-LAYER linear stability (bifurcation) analysis of a compressed
incompressible neo-Hookean laminate on a rigid base, with an independent
surface tension on every internal interface and on the free top surface.

This generalizes scripts/bilayer_stability.py (2 layers) to an arbitrary
stack. Layer k (k = 1..N, bottom to top) has shear modulus mu_k and
thickness H_k; it occupies Y in [Y_{k-1}, Y_k], Y_0 = 0 (rigid base),
Y_k = sum_{j<=k} H_j. Interface i (i = 1..N-1) sits at Y_i with tension
gamma_int[i-1]; the free top at Y_N carries gamma_top.

All layers share the common in-plane stretch lambda (bonded, uniform
compression), so within each layer the vertical envelope f2 obeys the
same 4th-order ODE and is a sum of four modes exp(+-K eta), exp(+-K eta
/lambda^2) in the layer-local coordinate eta. Each layer contributes 4
unknowns -> 4N total.

Conditions (4N):
  base (Y_0):           f1=0, f2=0                                  (2)
  each interface (Y_i): f1, f2, t1 continuous; t2 jump = b_i f2     (4 each)
  top (Y_N):            t1=0, t2 + b_top f2 = 0                     (2)

with b = gamma K^2 / lambda. det M(lambda,K) = 0 is the bifurcation
condition. The null vector at threshold gives the mode shape, which we
use to CLASSIFY the first instability as a surface mode or a buried
interface mode (the question: can an interface buckle before a
surface-tension-stiffened top?).

Tractions on a horizontal face (same reduction as the bilayer/paper):
  t1 = -mu[(lam^2/K) f2'' + (K/lam^2) f2]
  t2 =  mu[(2+lam^4) f2'  - (lam^4/K^2) f2''']
"""

import argparse
import numpy as np

# raw exponential modes overflow for large K*H/lam^2 (thick layers); those scan
# points become non-finite det and are skipped as non-roots. Suppress once,
# module-wide, instead of per-call (the per-call context manager dominated runtime).
np.seterr(over="ignore", invalid="ignore", divide="ignore")


def basis(K, lam, eta):
    r = np.array([K, -K, K / lam**2, -K / lam**2])
    e = np.exp(r * eta)
    return e, r * e, r**2 * e, r**3 * e   # f2, f2', f2'', f2'''


def fields(K, lam, mu, eta):
    f2, f2p, f2pp, f2ppp = basis(K, lam, eta)
    f1 = -(lam**2 / K) * f2p
    t1 = -mu * ((lam**2 / K) * f2pp + (K / lam**2) * f2)
    t2 = mu * ((2.0 + lam**4) * f2p - (lam**4 / K**2) * f2ppp)
    return f1, f2, t1, t2


def assemble(lam, K, mus, Hs, g_int, g_top, base="fixed"):
    """Build the 4N x 4N boundary/interface matrix.
    base="fixed":   bonded substrate, f1=0 and f2=0 (no slip).
    base="sliding": frictionless roller, t1=0 and f2=0 (free in-plane slip);
                    this matches the FE roller base used in the mechanical sweep."""
    N = len(mus)
    M = np.zeros((4 * N, 4 * N))
    row = 0

    def cols(k):  # column slice for layer k (0-based)
        return slice(4 * k, 4 * k + 4)

    # base: layer 0 at eta=0
    f1, f2, t1, t2 = fields(K, lam, mus[0], 0.0)
    M[row, cols(0)] = (t1 if base == "sliding" else f1); row += 1   # sliding: shear-free
    M[row, cols(0)] = f2; row += 1                                   # both: no normal displacement

    # interfaces
    for i in range(N - 1):
        f1a, f2a, t1a, t2a = fields(K, lam, mus[i], Hs[i])      # top of layer i
        f1b, f2b, t1b, t2b = fields(K, lam, mus[i + 1], 0.0)    # bottom of layer i+1
        b = g_int[i] * K**2 / lam
        M[row, cols(i)] = f1a;  M[row, cols(i + 1)] = -f1b;  row += 1   # f1 cont
        M[row, cols(i)] = f2a;  M[row, cols(i + 1)] = -f2b;  row += 1   # f2 cont
        M[row, cols(i)] = t1a;  M[row, cols(i + 1)] = -t1b;  row += 1   # t1 cont
        # t2 jump: t2_{i+1} - t2_i = b f2  -> attach b to upper (i+1) side
        M[row, cols(i)] = -t2a; M[row, cols(i + 1)] = t2b - b * f2b; row += 1

    # top: layer N-1 at eta=H_{N-1}
    f1, f2, t1, t2 = fields(K, lam, mus[-1], Hs[-1])
    bt = g_top * K**2 / lam
    M[row, cols(N - 1)] = t1; row += 1
    M[row, cols(N - 1)] = t2 + bt * f2; row += 1

    # row scaling for conditioning
    for r in range(4 * N):
        s = np.max(np.abs(M[r]))
        if s > 0:
            M[r] /= s
    return M


def detM(lam, K, mus, Hs, g_int, g_top, base="fixed"):
    return np.linalg.det(assemble(lam, K, mus, Hs, g_int, g_top, base))


def crit_lambda_at_K(mus, Hs, g_int, g_top, K, lam_hi=0.999, lam_lo=0.08, n=240, base="fixed"):
    lams = np.linspace(lam_hi, lam_lo, n)
    prev = detM(lams[0], K, mus, Hs, g_int, g_top, base)
    for j in range(1, len(lams)):
        cur = detM(lams[j], K, mus, Hs, g_int, g_top, base)
        if np.isfinite(prev) and np.isfinite(cur) and prev * cur < 0:
            a, b = lams[j], lams[j - 1]
            fa = cur
            for _ in range(40):
                m = 0.5 * (a + b)
                fm = detM(m, K, mus, Hs, g_int, g_top, base)
                if fa * fm < 0:
                    b = m
                else:
                    a, fa = m, fm
            return 0.5 * (a + b)
        prev = cur
    return None


def critical(mus, Hs, g_int, g_top, Ks, refine=True, base="fixed", lam_lo=0.08):
    """Critical strain minimized over K. With refine=True, parabolically
    interpolate K_c from the discrete minimum's neighbours so the selected
    wavelength is smooth (not quantized to the K grid).  base selects the
    substrate condition ('fixed' bonded, or 'sliding' roller)."""
    eps = np.full(len(Ks), np.inf)
    lams = np.full(len(Ks), np.nan)
    for j, K in enumerate(Ks):
        lam = crit_lambda_at_K(mus, Hs, g_int, g_top, K, lam_lo=lam_lo, base=base)
        if lam is not None:
            eps[j] = 1.0 - lam; lams[j] = lam
    i = int(np.argmin(eps))
    if not np.isfinite(eps[i]):
        return None
    best = (eps[i], Ks[i], lams[i])
    if refine and 0 < i < len(Ks) - 1 and np.isfinite(eps[i - 1]) and np.isfinite(eps[i + 1]):
        e0, e1, e2 = eps[i - 1], eps[i], eps[i + 1]
        denom = (e0 - 2 * e1 + e2)
        if denom > 0:
            # parabola vertex in LOG-K: the spacing is uniform for a geomspace grid,
            # so the vertex formula is exact and K_c(gbar) comes out smooth (a linear-K
            # vertex formula on a log grid mis-locates K_c and makes l_c jagged).
            u0, u1, u2 = np.log(Ks[i - 1]), np.log(Ks[i]), np.log(Ks[i + 1])
            uv = u1 + 0.5 * (e0 - e2) / denom * (u2 - u1)
            kv = float(np.exp(uv))
            lamv = crit_lambda_at_K(mus, Hs, g_int, g_top, kv, lam_lo=lam_lo, base=base)
            if lamv is not None and (1.0 - lamv) <= best[0]:
                best = (1.0 - lamv, kv, lamv)
    return best


def mode_shape(lam, K, mus, Hs, g_int, g_top, npts=20):
    """Null vector -> vertical-displacement envelope f2(Y) sampled per layer.
    Returns (Y_array, f2_array) and the interface/surface deflection magnitudes."""
    M = assemble(lam, K, mus, Hs, g_int, g_top)
    _, s, vh = np.linalg.svd(M)
    c = vh[-1]                       # null-space vector (smallest singular value)
    Ys, F2 = [], []
    Y0 = 0.0
    bound_defl = []                  # |f2| at base, each interface, top
    for k in range(len(mus)):
        etas = np.linspace(0, Hs[k], npts)
        ck = c[4 * k:4 * k + 4]
        for eta in etas:
            f2 = basis(K, lam, eta)[0] @ ck
            Ys.append(Y0 + eta); F2.append(f2)
        # boundary deflection at top of this layer
        f2_top = basis(K, lam, Hs[k])[0] @ ck
        bound_defl.append(abs(f2_top))
        Y0 += Hs[k]
    F2 = np.array(F2)
    nrm = np.max(np.abs(F2)) or 1.0
    return np.array(Ys), F2 / nrm, np.array(bound_defl) / nrm


def parse_stack(spec):
    """spec like 'mu:H,mu:H,...' bottom->top."""
    mus, Hs = [], []
    for tok in spec.split(","):
        m, h = tok.split(":")
        mus.append(float(m)); Hs.append(float(h))
    return mus, Hs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stack", default="1:2,5:2",
                    help="'mu:H,mu:H,...' bottom to top (default bilayer 1:2,5:2)")
    ap.add_argument("--g_int", default=None,
                    help="comma list of interface tensions (len N-1); default all 0")
    ap.add_argument("--g_top", type=float, default=0.0)
    ap.add_argument("--Kmax", type=float, default=4.0)
    ap.add_argument("--nK", type=int, default=240)
    ap.add_argument("--mode", action="store_true", help="print mode shape classification")
    args = ap.parse_args()

    mus, Hs = parse_stack(args.stack)
    N = len(mus)
    if args.g_int is None:
        g_int = [0.0] * (N - 1)
    else:
        g_int = [float(x) for x in args.g_int.split(",")]
        assert len(g_int) == N - 1, f"need {N-1} interface tensions"

    Ks = np.linspace(0.05, args.Kmax, args.nK)
    res = critical(mus, Hs, g_int, g_top=args.g_top, Ks=Ks)
    Htot = sum(Hs)
    print(f"N={N} stack mu={mus} H={Hs}  g_int={g_int}  g_top={args.g_top}")
    if not res:
        print("  no bifurcation found in scan range"); return
    eps, K, lam = res
    print(f"  eps_c={eps:.4f}  lambda_c={lam:.4f}  K_c={K:.3f}  l_c={2*np.pi/K:.2f}  l_c/H={2*np.pi/K/Htot:.2f}")

    if args.mode:
        Y, f2, bd = mode_shape(lam, K, mus, Hs, g_int, args.g_top)
        names = ["interface %d (Y=%.2f)" % (i + 1, sum(Hs[:i + 1])) for i in range(N - 1)]
        names.append("TOP surface (Y=%.2f)" % Htot)
        kmax = int(np.argmax(bd))
        print("  boundary deflection magnitudes (normalized):")
        for nm, d in zip(names, bd):
            star = "  <== first instability localizes here" if names.index(nm) == kmax else ""
            print(f"     {nm:28} |f2|={d:.3f}{star}")


if __name__ == "__main__":
    main()
