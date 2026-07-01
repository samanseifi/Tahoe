#!/usr/bin/env python3
r"""
Linear stability (bifurcation) analysis of a compressed incompressible
neo-Hookean BILAYER on a rigid base, with surface tension on the free top
surface (gamma_top) and on the internal interface (gamma_int).

This is the multilayer generalization of the single-layer perturbation
analysis in paper_draft/main.tex (eqs. ODE_2 and matrix A).

Geometry (reference):
    rigid base ----- Y = 0
      layer 1 (substrate, mu1)   0   <= Y <= H1
    interface ------ Y = H1   (gamma_int)
      layer 2 (film, mu2)        H1  <= Y <= H1+H2
    free top ------- Y = H1+H2   (gamma_top)

Base state: uniform plane-strain incompressible compression, in-plane
stretch lambda (epsilon = 1 - lambda), F0 = diag(lambda, 1/lambda) in both
layers (bonded, common in-plane stretch). Per-layer base pressure
pi0_k = mu_k / lambda^2 (traction-free top, flat interface).

Perturbation (reference coords X, Y), wavenumber K = 2*pi/L:
    dx1 = f1(Y) sin(K X),  dx2 = f2(Y) cos(K X),  dpi = f3(Y) cos(K X)

In each layer f2 satisfies   f2'''' - K^2 (1+lambda^-4) f2'' + K^4 lambda^-4 f2 = 0,
so f2 = sum_m c_m exp(r_m Y),  r_m in {+K, -K, +K/lambda^2, -K/lambda^2}.

Incompressibility  ->  f1 = -(lambda^2/K) f2'
Pressure (from horizontal equilibrium) -> f3 = mu lambda^3 (f2'''/K^2 - f2')
Nominal tractions on a horizontal face (normal = e_Y):
    t1 (shear)  = -mu [ (lambda^2/K) f2'' + (K/lambda^2) f2 ]
    t2 (normal) =  mu [ (2 + lambda^4) f2' - (lambda^4/K^2) f2''' ]

Surface tension enters the NORMAL traction BC as a restoring Laplace
pressure. Current surface wavenumber k = K/lambda; Laplace pressure
increment ~ gamma k^2 w; pushed to the reference face (Nanson factor
|F0^-T e_Y| = lambda) gives a restoring nominal traction
    t2_surf = (gamma K^2 / lambda) f2      (stabilizing).

8 conditions / 8 unknowns (c1[0:4], c2[0:4]):
  base   : f1_1(0) = 0,  f2_1(0) = 0
  interf : f1 cont, f2 cont, t1 cont, [t2] = (gamma_int K^2/lambda) f2
  top    : t1_2(top) = 0,  t2_2(top) + (gamma_top K^2/lambda) f2_2(top) = 0

det M(lambda, K; ...) = 0 is the bifurcation condition. For each K we find
the largest lambda (smallest strain) with a sign change; epsilon_c = min
over K of (1 - lambda), and K_c is the selected wavenumber.

Validation (run with --validate): single layer (mu2=mu1, gamma=0) must
recover Biot epsilon ~ 0.456 (lambda ~ 0.544) in the thick-layer limit.
"""

import argparse
import numpy as np


def layer_basis(K, lam, eta):
    """Return f2 and derivatives' coefficient row-vectors at local coord eta.
    Basis exponents r = [K, -K, K/lam^2, -K/lam^2]; returns dict of
    length-4 arrays giving [f2, f2', f2'', f2'''] contributions per mode."""
    r = np.array([K, -K, K / lam**2, -K / lam**2])
    e = np.exp(r * eta)
    return {
        0: e,            # f2
        1: r * e,        # f2'
        2: r**2 * e,     # f2''
        3: r**3 * e,     # f2'''
    }


def tractions(K, lam, mu, d):
    """Given derivative rows d[0..3] (each length-4), return row-vectors for
    f1, f2, t1 (shear), t2 (normal)."""
    f2, f2p, f2pp, f2ppp = d[0], d[1], d[2], d[3]
    f1 = -(lam**2 / K) * f2p
    t1 = -mu * ((lam**2 / K) * f2pp + (K / lam**2) * f2)
    t2 = mu * ((2.0 + lam**4) * f2p - (lam**4 / K**2) * f2ppp)
    return f1, f2, t1, t2


def detM_bilayer(lam, K, mu1, mu2, H1, H2, g_top, g_int):
    """Assemble the 8x8 boundary/interface matrix and return its determinant.
    Unknowns: c1 (layer1, local eta in [0,H1]), c2 (layer2, local eta in [0,H2])."""
    M = np.zeros((8, 8))
    b = g_int * K**2 / lam   # interfacial surface-tension coefficient
    t = g_top * K**2 / lam   # top surface-tension coefficient

    # layer 1 at base (eta=0) and interface (eta=H1)
    d1_0 = layer_basis(K, lam, 0.0)
    d1_H = layer_basis(K, lam, H1)
    f1_1_0, f2_1_0, t1_1_0, t2_1_0 = tractions(K, lam, mu1, d1_0)
    f1_1_H, f2_1_H, t1_1_H, t2_1_H = tractions(K, lam, mu1, d1_H)

    # layer 2 at interface (eta=0) and top (eta=H2)
    d2_0 = layer_basis(K, lam, 0.0)
    d2_T = layer_basis(K, lam, H2)
    f1_2_0, f2_2_0, t1_2_0, t2_2_0 = tractions(K, lam, mu2, d2_0)
    f1_2_T, f2_2_T, t1_2_T, t2_2_T = tractions(K, lam, mu2, d2_T)

    # rows (cols 0:4 -> c1, cols 4:8 -> c2)
    M[0, 0:4] = f1_1_0                       # base: f1=0
    M[1, 0:4] = f2_1_0                       # base: f2=0
    M[2, 0:4] = f1_1_H;  M[2, 4:8] = -f1_2_0  # interface: f1 continuous
    M[3, 0:4] = f2_1_H;  M[3, 4:8] = -f2_2_0  # interface: f2 continuous
    M[4, 0:4] = t1_1_H;  M[4, 4:8] = -t1_2_0  # interface: shear traction continuous
    # interface: t2^+ - t2^- = b f2  (restoring Laplace tension).
    # Vacuum-above limit (t2^+=0) -> t2^- + b f2 = 0, identical to the
    # validated top-surface condition. b attached to layer-2 (upper) side.
    M[5, 0:4] = -t2_1_H; M[5, 4:8] = t2_2_0 - b * f2_2_0
    M[6, 4:8] = t1_2_T                       # top: shear-free
    M[7, 4:8] = t2_2_T + t * f2_2_T          # top: surface tension

    # scale rows to keep determinant well-conditioned
    for i in range(8):
        s = np.max(np.abs(M[i]))
        if s > 0:
            M[i] /= s
    return np.linalg.det(M)


def detM_single(lam, K, mu, H, g_top):
    """Single homogeneous layer on rigid base, free top + surface tension.
    4 unknowns. Used for validation against Biot."""
    M = np.zeros((4, 4))
    t = g_top * K**2 / lam
    d0 = layer_basis(K, lam, 0.0)
    dH = layer_basis(K, lam, H)
    f1_0, f2_0, _, _ = tractions(K, lam, mu, d0)
    _, f2_H, t1_H, t2_H = tractions(K, lam, mu, dH)
    M[0] = f1_0
    M[1] = f2_0
    M[2] = t1_H
    M[3] = t2_H + t * f2_H
    for i in range(4):
        s = np.max(np.abs(M[i]))
        if s > 0:
            M[i] /= s
    return np.linalg.det(M)


def critical_lambda_at_K(detfun, K, lam_hi=0.999, lam_lo=0.30, n=400):
    """Scan lambda downward; return the largest lambda (smallest strain) where
    det changes sign (first bifurcation), or None."""
    lams = np.linspace(lam_hi, lam_lo, n)
    prev = detfun(lams[0], K)
    for lam in lams[1:]:
        cur = detfun(lam, K)
        if np.isfinite(prev) and np.isfinite(cur) and prev * cur < 0:
            # bisect
            a, b = lam, lam + (lams[0] - lams[1])
            fa = cur
            for _ in range(60):
                m = 0.5 * (a + b)
                fm = detfun(m, K)
                if fa * fm < 0:
                    b = m
                else:
                    a, fa = m, fm
            return 0.5 * (a + b)
        prev = cur
    return None


def critical_strain(detfun, Ks):
    """Minimize critical strain over wavenumber set Ks."""
    best = None
    for K in Ks:
        lam = critical_lambda_at_K(detfun, K)
        if lam is None:
            continue
        eps = 1.0 - lam
        if best is None or eps < best[0]:
            best = (eps, K, lam)
    return best  # (eps_c, K_c, lam_c) or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true",
                    help="recover Biot single-layer limit")
    ap.add_argument("--H1", type=float, default=2.0)
    ap.add_argument("--H2", type=float, default=2.0)
    ap.add_argument("--mu1", type=float, default=1.0)
    ap.add_argument("--mu2", type=float, default=5.0)
    ap.add_argument("--g_top", type=float, default=0.0)
    ap.add_argument("--g_int", type=float, default=0.0)
    args = ap.parse_args()

    if args.validate:
        # thick single layer, gamma=0 -> Biot lambda~0.544 (eps~0.456)
        H = 50.0
        Ks = np.linspace(0.02, 2.0, 120) / 1.0
        det = lambda lam, K: detM_single(lam, K, 1.0, H, 0.0)
        res = critical_strain(det, Ks)
        print("== Validation: single layer, gamma=0, thick (Biot target eps=0.456) ==")
        if res:
            print(f"   eps_c = {res[0]:.4f}  (lambda_c = {res[2]:.4f}),  K_c*H = {res[1]*H:.2f}")
        else:
            print("   no bifurcation found")
        # surface-tension shift check (paper: eps_c rises with gamma_bar)
        for gb in (0.0, 0.5, 2.0):
            Hf = 4.0
            g = gb * 1.0 * Hf  # gamma = gbar*mu*H
            Ks2 = np.linspace(0.05, 3.0, 160) / 1.0
            det2 = lambda lam, K: detM_single(lam, K, 1.0, Hf, g)
            r = critical_strain(det2, Ks2)
            if r:
                print(f"   gamma_bar={gb:>4}: eps_c={r[0]:.3f}, l_c/H={2*np.pi/r[1]/Hf:.2f}")
        return

    print(f"== Bilayer: mu2/mu1={args.mu2/args.mu1:g}, H1={args.H1}, H2={args.H2} ==")
    Ks = np.linspace(0.05, 4.0, 200) / 1.0
    for label, g_int in (("interface OFF", 0.0), (f"interface gamma={args.g_int}", args.g_int)):
        gi = g_int
        det = lambda lam, K: detM_bilayer(lam, K, args.mu1, args.mu2,
                                          args.H1, args.H2, args.g_top, gi)
        res = critical_strain(det, Ks)
        if res:
            eps, K, lam = res
            print(f"   {label:>28}: eps_c={eps:.4f}, lambda_c={lam:.4f}, "
                  f"K_c={K:.3f}, l_c={2*np.pi/K:.2f}")
        else:
            print(f"   {label:>28}: no bifurcation found in scan range")


if __name__ == "__main__":
    main()
