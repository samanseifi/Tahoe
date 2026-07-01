#!/usr/bin/env python3
r"""
Weakly nonlinear Landau coefficient for surface wrinkling of a compressed
incompressible neo-Hookean film WITH surface tension (gam>0, where the mode is
isolated -- the gam=0 half-space is scale-free/degenerate, see weakly_nonlinear.py).

Nonlinear forcing is computed by NUMERICAL FFT PROJECTION of the exact residual
(bulk is exactly quadratic: S2 = -p cof(H); constraint cof(F0):H + detH = 0),
which avoids hand harmonic algebra.  Amplitude A = surface vertical displacement
(f2(0) normalized to 1).  Result:  lam = lam_c + LAM2 * A^2, i.e.
  eps - eps_c = -LAM2 * A^2.    LAM2 > 0 -> subcritical (crease-like).

Stage A: bulk geometric nonlinearity + linear Young-Laplace (in operator).
Stage B (add_surface_nl=True): + nonlinear YL curvature on the top BC.
"""
import numpy as np
import weakly_nonlinear as wn


def _cof(H):
    # H[i][j] each (n,M);  cof = [[H11,-H10],[-H01,H00]]
    return [[H[1][1], -H[1][0]], [-H[0][1], H[0][0]]]


def landau(N, H, lam, k0, gam, M=24, add_surface_nl=False):
    X2, D, D2 = wn.operators(N, H)
    n = N + 1
    X1 = (2 * np.pi / k0) * np.arange(M) / M
    kk = np.fft.fftfreq(M, d=1.0 / M) * k0          # X1 wavenumbers

    def dX1(F):
        return np.real(np.fft.ifft(1j * kk[None, :] * np.fft.fft(F, axis=1), axis=1))

    def dX2(F):
        return D @ F

    def proj(F, m):
        """cos/sin envelopes of harmonic m (length n)."""
        Fh = np.fft.fft(F, axis=1)
        if m == 0:
            return np.real(Fh[:, 0]) / M, np.zeros(n)
        return 2 * np.real(Fh[:, m]) / M, -2 * np.imag(Fh[:, m]) / M

    def make(env_c, env_s, m):
        cm = np.cos(m * k0 * X1); sm = np.sin(m * k0 * X1)
        return np.outer(env_c, cm) + np.outer(env_s, sm)

    def grads(wx, wy):
        return [[dX1(wx), dX2(wx)], [dX1(wy), dX2(wy)]]

    def S2(p, Hg):
        cf = _cof(Hg)
        return [[-p * cf[i][j] for j in (0, 1)] for i in (0, 1)]   # S2[i][j]

    def divS(St):
        return [dX1(St[i][0]) + dX2(St[i][1]) for i in (0, 1)]

    def detH(Hg):
        return Hg[0][0] * Hg[1][1] - Hg[0][1] * Hg[1][0]

    # ---------- order 1 ----------
    f1, f2, f3 = wn.eigenmode(N, D, D2, lam, k0, gam)
    s0 = f2[0]; f1, f2, f3 = f1 / s0, f2 / s0, f3 / s0
    w1x = make(np.zeros(n), f1, 1)
    w1y = make(f2, np.zeros(n), 1)
    p1 = make(f3, np.zeros(n), 1)
    H1 = grads(w1x, w1y)
    S2_11 = S2(p1, H1)
    dS2_1 = divS(S2_11)                                # (DivS2)_i for w1*w1

    # ---------- order 2 forcing: L w2 = -DivS2 ; constraint cof(F0):H2 = -detH1
    detH1 = detH(H1)

    def solve_order2(m):
        rI_c, _ = proj(-dS2_1[0], m)   # comp1 eq carried by sin; use cos? handle parity below
        rI_s, _ = (None, None)
        # (DivS2)_1 is sin-type (odd) -> use sin envelope; (DivS2)_2 cos-type (even)
        cI, sI = proj(dS2_1[0], m)     # comp1
        cII, sII = proj(dS2_1[1], m)   # comp2
        cC, sC = proj(detH1, m)        # constraint
        # operator rows: EI carried by sin (comp1), EII by cos (comp2), EIII by cos
        rhsI = -sI                     # = -(DivS2)_1   (sin envelope)
        rhsII = -cII                   # = -(DivS2)_2   (cos envelope)
        rhsIII = -cC                   # = -detH1       (cos)
        # surface traction forcing  S1[i,1](w2) = -S2[i,1](w1)  at node 0
        cS01, sS01 = proj(S2_11[0][1], m)   # shear S2[0,1]
        cS11, sS11 = proj(S2_11[1][1], m)   # normal S2[1,1]
        bc_shear = -sS01[0]            # shear carried by sin
        bc_norm = -cS11[0]             # normal carried by cos
        g1, g2, g3 = wn._forced_solve(N, D, D2, lam, m * k0, gam,
                                      rhsI, rhsII, rhsIII, bc_shear, bc_norm)
        return g1, g2, g3

    # mean (m=0) and second harmonic (m=2)
    h1, h2, h3 = solve_order2(0)
    g1, g2, g3 = solve_order2(2)
    w2x = make(np.zeros(n), h1, 0) + make(g1 * 0, g1, 2)   # comp1: mean(odd->0 at k=0)+2nd sin
    # NOTE mean comp1 should be ~0; keep general:
    w2x = np.outer(h1, np.ones(M)) + make(np.zeros(n), g1, 2)
    w2y = np.outer(h2, np.ones(M)) + make(g2, np.zeros(n), 2)
    p2 = np.outer(h3, np.ones(M)) + make(g3, np.zeros(n), 2)
    H2 = grads(w2x, w2y)

    # ---------- order 3 resonant forcing (k0 harmonic) ----------
    # bulk: A^3 part of -DivS2(w) = -Div(p1 cofH2 + p2 cofH1)
    cross = [[-(p1 * _cof(H2)[i][j] + p2 * _cof(H1)[i][j]) for j in (0, 1)] for i in (0, 1)]
    dCross = divS(cross)
    # constraint A^3: cof(F0):H3 = -(cofH1:H2)   where cofH1:H2 = sum cof(H1)[i][j] H2[i][j]
    cf1 = _cof(H1)
    cofH1_H2 = sum(cf1[i][j] * H2[i][j] for i in (0, 1) for j in (0, 1))
    # project to k0 (m=1): comp1 sin, comp2 cos, constraint cos
    _, sI3 = proj(dCross[0], 1)
    cII3, _ = proj(dCross[1], 1)
    cC3, _ = proj(cofH1_H2, 1)
    rhsI3 = -sI3
    rhsII3 = -cII3
    rhsIII3 = -cC3
    # surface order-3 traction: S1[i,1](w3) = -(cross)[i,1] at node0 (+ YL cubic if Stage B)
    _, sCr01 = proj(cross[0][1], 1)
    cCr11, _ = proj(cross[1][1], 1)
    bc_shear3 = -sCr01[0]
    bc_norm3 = -cCr11[0]
    if add_surface_nl:
        # nonlinear Young-Laplace: kappa = -eta'' (1 - 3/2 eta'^2); cubic part (3/2)eta1'' eta1'^2
        # surface height (current) eta(x1)= w1y(0)*A ; x1 = lam X1 + w1x(0)*A
        # leading cubic curvature correction -> nominal normal traction forcing at k0.
        eta1 = f2[0] * np.cos(k0 * X1)          # = cos(k0 X1) (f2[0]=1)
        etap = -k0 * np.sin(k0 * X1)            # d eta1/d x1approx /lam
        etapp = -k0**2 * np.cos(k0 * X1)
        cub = 1.5 * etapp * (etap / lam)**2     # (3/2) eta'' eta'^2 (current coords ~ /lam)
        cc, _ = proj(np.outer(np.ones(n), cub), 1)
        bc_norm3 += (gam / lam) * cc[0]         # gamma * cubic curvature -> nominal (/lam area)

    # ---------- solvability ----------
    adj = wn.adjoint_mode(N, D, D2, lam, k0, gam)
    b3 = np.zeros(3 * n)
    b3[0:n] = rhsI3; b3[n:2 * n] = rhsII3; b3[2 * n:3 * n] = rhsIII3
    b3[0] = bc_shear3; b3[n] = bc_norm3
    num = adj @ b3
    e1 = np.concatenate([f1, f2, f3])
    dl = 1e-3
    Lp = wn.linear_matrix(N, D, D2, lam + dl, k0, gam)
    Lm = wn.linear_matrix(N, D, D2, lam - dl, k0, gam)
    den = adj @ (((Lp - Lm) / (2 * dl)) @ e1)
    LAM2 = -num / den
    # second-order field magnitude (sanity: should be O(1), not blowing up)
    mag2 = max(np.abs(g2).max(), np.abs(h2).max())
    return LAM2, mag2, num, den


if __name__ == "__main__":
    N = 48; H = 4.0
    print("Weakly nonlinear Landau coefficient (eps - eps_c = -LAM2 A^2):")
    print("  LAM2>0 => SUBCRITICAL (crease-like) ; LAM2<0 => SUPERCRITICAL (smooth wrinkle)")
    X2, D, D2 = wn.operators(N, H)
    for gam in (1.0, 2.0, 4.0, 8.0):
        # find selected k0 and lam_c
        best = None
        for k in np.linspace(0.3, 2.2, 26):
            lc, sv = wn.critical_lambda(N, D, D2, k, gam)
            if lc is not None and (best is None or (1 - lc) < best[0]):
                best = (1 - lc, k, lc)
        if best is None:
            print(f"  gbar={gam/H:.2f}: no isolated mode found"); continue
        eps, k0, lc = best
        L2, mag2, num, den = landau(N, H, lc, k0, gam)
        print(f"  gbar={gam/H:.2f} (gam={gam}): k0={k0:.2f} eps_c={eps:.3f} | "
              f"LAM2={L2:+.3f}  (|w2|~{mag2:.1f})  [{'SUB' if L2>0 else 'SUPER'}critical]", flush=True)
