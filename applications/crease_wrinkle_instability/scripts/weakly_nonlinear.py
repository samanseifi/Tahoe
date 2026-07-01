#!/usr/bin/env python3
r"""
Weakly nonlinear (Landau) analysis of surface wrinkling of a compressed
incompressible neo-Hookean film, built on the SAME incremental theory as the
linear analysis (elec_stability.py) but carried to third order in the modal
amplitude A.  Mechanical case first (E=0); the electric coupling plugs into the
top boundary condition afterwards.

Reference-config exact kinematics (2D, plane strain), F = F0 + H, F0=diag(lam,1/lam):
  S   = mu F - pi cof(F),    cof(F) = cof(F0) + cof(H)   (EXACTLY linear in H)
  det F = 1  ->  cof(F0):H + det H = 0                   (constraint, exactly quadratic)
=> bulk nonlinearity is PURELY QUADRATIC:
  S1 = mu H - pi0 cof(H) - p cof(F0)      (linear)
  S2 = - p cof(H)                          (quadratic; the only nonlinear stress)
Equilibrium Div S = 0 order by order; pi0 = mu/lam^2.

Per Fourier harmonic in X1 the incremental fields reduce (mu=1) to ODEs in X2:
  w1 = f1(X2) sin(kX1), w2 = f2(X2) cos(kX1), p = f3(X2) cos(kX1)
  (I)   f1'' - k^2 f1 + (k/lam) f3 = 0
  (II)  f2'' - k^2 f2 - lam f3'    = 0
  (III) (k/lam) f1 + lam f2'       = 0      (incompressibility)
Top surface X2=0 traction (outward normal -e2), free + Young-Laplace:
  shear : f1'(0) - pi0 k f2(0) = 0
  normal: f2'(0) - pi0 k f1(0) - lam f3(0) + (gam k^2/lam) f2(0) = 0
Base X2=H clamped: f1(H)=f2(H)=0.

This module (stage 1) sets up the LINEAR solver and validates Biot (eps_c=0.456,
gam=0, large kH) and the surface-tension shift.  Stages 2-3 (second harmonic +
mean field, third-order solvability -> Landau coefficient) are added next.
"""
import numpy as np


def cheb(N):
    """Chebyshev differentiation matrix on Gauss-Lobatto nodes x in [1,-1]."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = (np.outer(c, 1.0 / c)) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def operators(N, H):
    """Chebyshev D, D2 on physical X2 in [0,H] (X2=0 at surface, X2=H base)."""
    D1, xc = cheb(N)              # xc in [1,-1]
    # map xc in [1,-1] -> X2 in [0,H]:  X2 = H*(1-xc)/2
    X2 = H * (1.0 - xc) / 2.0
    fac = -2.0 / H               # dxc/dX2
    D = fac * D1
    D2 = D @ D
    return X2, D, D2


def linear_matrix(N, D, D2, lam, k, gam):
    """Assemble 3(N+1) system for (f1,f2,f3); rows replaced by BCs.  mu=1."""
    n = N + 1
    pi0 = 1.0 / lam**2
    I = np.eye(n)
    Z = np.zeros((n, n))
    # interior equations (we will overwrite boundary rows with BCs)
    # (I)  f1'' - k^2 f1 + (k/lam) f3 = 0
    EI = np.hstack([D2 - k**2 * I, Z, (k / lam) * I])
    # (II) f2'' - k^2 f2 - lam f3'    = 0
    EII = np.hstack([Z, D2 - k**2 * I, -lam * D])
    # (III)(k/lam) f1 + lam f2'       = 0
    EIII = np.hstack([(k / lam) * I, lam * D, Z])
    M = np.vstack([EI, EII, EIII])
    # node indexing: block f1 = rows[0:n], f2=[n:2n], f3=[2n:3n]
    # surface is node index 0 (X2=0), base is node index N (X2=H)
    s, b = 0, N
    # --- replace rows with BCs ---
    # use equation-(I) row at surface (s) and base (b) for f1 BCs; etc. Simpler:
    # overwrite specific equation rows that touch the boundary nodes.
    def row(blk, node):
        return blk * n + node
    # Top shear: f1'(0) - pi0 k f2(0) = 0   -> put in EI surface row
    M[row(0, s), :] = np.hstack([D[s, :], -pi0 * k * I[s, :], Z[s, :]])
    # Top normal: f2'(0) - pi0 k f1(0) - lam f3(0) - (gam k^2/lam) f2(0) = 0 -> EII surface row
    # (the Young-Laplace term is stabilizing; sign fixed by matching eps_c rising with gamma)
    M[row(1, s), :] = np.hstack([-pi0 * k * I[s, :],
                                 D[s, :] - (gam * k**2 / lam) * I[s, :],
                                 -lam * I[s, :]])
    # Base clamp f1(H)=0, f2(H)=0 -> EI base row, EII base row
    M[row(0, b), :] = np.hstack([I[b, :], Z[b, :], Z[b, :]])
    M[row(1, b), :] = np.hstack([Z[b, :], I[b, :], Z[b, :]])
    return M


def smallest_sv(N, D, D2, lam, k, gam):
    M = linear_matrix(N, D, D2, lam, k, gam)
    return np.linalg.svd(M, compute_uv=False)[-1]


def critical_lambda(N, D, D2, k, gam, lo=0.3, hi=0.98, ncoarse=160, tol=2e-3):
    """largest lam (least compression) at a genuine interior singular dip of sv(lam).
    Returns (lam_c, sv) or (None, None) if no true neutral mode in (lo,hi)."""
    lams = np.linspace(hi, lo, ncoarse)
    sv = np.array([smallest_sv(N, D, D2, l, k, gam) for l in lams])
    # interior local minima (dips) below tolerance = genuine neutral modes
    cand = [i for i in range(1, len(lams) - 1)
            if sv[i] < sv[i - 1] and sv[i] < sv[i + 1] and sv[i] < tol]
    if not cand:
        return None, None
    i = cand[0]                      # largest lam (lams descending) = first instability
    a, b = lams[i + 1], lams[i - 1]
    for _ in range(60):
        m1 = a + (b - a) / 3; m2 = b - (b - a) / 3
        if smallest_sv(N, D, D2, m1, k, gam) < smallest_sv(N, D, D2, m2, k, gam):
            b = m2
        else:
            a = m1
    lc = 0.5 * (a + b)
    return lc, smallest_sv(N, D, D2, lc, k, gam)


def eigenmode(N, D, D2, lam, k, gam):
    """right null vector (f1,f2,f3) of the linear operator at a neutral point."""
    M = linear_matrix(N, D, D2, lam, k, gam)
    U, S, Vt = np.linalg.svd(M)
    v = Vt[-1]
    n = N + 1
    return v[:n], v[n:2 * n], v[2 * n:]      # f1,f2,f3 envelopes


def adjoint_mode(N, D, D2, lam, k, gam):
    """left null vector (adjoint) of the linear operator at a neutral point."""
    M = linear_matrix(N, D, D2, lam, k, gam)
    U, S, Vt = np.linalg.svd(M)
    return U[:, -1]                           # 3n-vector


if __name__ == "__main__":
    N = 40
    print("== LINEAR validation ==", flush=True)
    # Biot half-space limit: gam=0, large kH -> eps_c -> 0.456 (use H=2, large k)
    for H, k in [(2.0, 8.0), (2.0, 12.0), (2.0, 16.0)]:
        X2, D, D2 = operators(N, H)
        lc, sv = critical_lambda(N, D, D2, k, 0.0)
        print(f"  gam=0  kH={k*H:.0f}:  lam_c={lc:.4f}  eps_c={1-lc:.4f}  (Biot 0.4563)  sv={sv:.1e}", flush=True)
    print("  -- surface-tension shift (finite film H=4, min over k) --", flush=True)
    H = 4.0; X2, D, D2 = operators(N, H)
    for gam in (0.0, 2.0, 4.0, 8.0):
        ks = np.linspace(0.4, 3.0, 20)
        best = None
        for k in ks:
            lc, sv = critical_lambda(N, D, D2, k, gam)
            if lc is not None and (best is None or (1 - lc) < best[0]):
                best = (1 - lc, k, lc)
        if best:
            print(f"  gbar={gam/(1.0*H):.2f} (gam={gam}): eps_c={best[0]:.3f} at k={best[1]:.2f} (l/H={2*np.pi/best[1]/H:.2f})", flush=True)


# =====================================================================
# Stage 2-3: weakly nonlinear Landau coefficient (mechanical, gam=0 first)
# Bulk is exactly quadratic: S2 = -p cof(H); constraint cof(F0):H + detH = 0.
# Harmonics: w1 ~ k0 ; second order -> mean(k=0) + 2nd harmonic(k=2k0).
# =====================================================================

def _forced_solve(N, D, D2, lam, k, gam, rhsI, rhsII, rhsIII, bc_shear, bc_norm):
    """Solve the collocation operator at wavenumber k with interior forcing
    (rhsI,rhsII,rhsIII) and top-surface traction forcing (bc_shear,bc_norm).
    Base clamped (0). Returns envelopes g1,g2,g3."""
    n = N + 1
    M = linear_matrix(N, D, D2, lam, k, gam)
    b = np.zeros(3 * n)
    b[0:n] = rhsI                       # eq I rows
    b[n:2 * n] = rhsII                  # eq II rows
    b[2 * n:3 * n] = rhsIII             # eq III rows
    s, bb = 0, N
    # overwrite BC rows' RHS (rows were replaced in linear_matrix at these indices)
    b[0 * n + s] = bc_shear             # top shear row (EI surface)
    b[1 * n + s] = bc_norm              # top normal row (EII surface)
    b[0 * n + bb] = 0.0                 # base f1(H)=0
    b[1 * n + bb] = 0.0                 # base f2(H)=0
    x = np.linalg.solve(M, b)
    return x[:n], x[n:2 * n], x[2 * n:]


def landau_gam0(N, H, lam, k0, dlam=1e-3):
    """Landau coefficient lam2 in  eps - eps_c = -lam2 A^2  for gam=0.
    Sign: lam2>0 => states at SMALLER compression (subcritical, crease-like)."""
    X2, D, D2 = operators(N, H)
    n = N + 1
    # --- order 1: eigenmode (normalize surface vertical displacement to 1) ---
    f1, f2, f3 = eigenmode(N, D, D2, lam, k0, 0.0)
    sc = f2[0]
    f1, f2, f3 = f1 / sc, f2 / sc, f3 / sc
    d = lambda y: D @ y                 # X2-derivative
    f1p, f2p, f3p = d(f1), d(f2), d(f3)

    # --- order 2 forcing (mean k=0 and 2nd harmonic k=2k0) ---
    # 2nd harmonic interior RHS  (operator = -(DivS2)_i ; constraint = -detH1)
    rI2  = -(k0 * f3 * f2p - 0.5 * k0 * d(f2 * f3))
    rII2 = -(k0 * f3 * f1p - 0.5 * k0 * d(f1 * f3))
    rIII2 = -0.5 * k0 * (f1 * f2p - f1p * f2)
    bcs2 = 0.5 * k0 * f2[0] * f3[0]     # shear forcing at surface
    bcn2 = 0.5 * k0 * f1[0] * f3[0]     # normal forcing at surface
    g1, g2, g3 = _forced_solve(N, D, D2, lam, 2 * k0, 0.0, rI2, rII2, rIII2, bcs2, bcn2)
    g1p, g2p, g3p = d(g1), d(g2), d(g3)

    # mean field k=0: only vertical w20_2 and p20 ; w20_1=0
    # constraint: lam*w20_2' = -0.5 k0 (f1 f2)'  ; equilibrium fixes p20.
    rI0  = np.zeros(n)
    rII0 = -0.5 * k0 * d(f1 * f3)        # (DivS2)_2 mean = -0.5k0 (f1 f3)'  -> operator = +...
    rIII0 = -0.5 * k0 * d(f1 * f2)
    h1, h2, h3 = _forced_solve(N, D, D2, lam, 0.0, 0.0, rI0, rII0, rIII0, 0.0, 0.0)
    h1p, h2p, h3p = d(h1), d(h2), d(h3)

    # --- order 3 resonant (k=k0) forcing from S3 = -(p1 cof(H2)+p2 cof(H1)) ---
    # Build third-order interior forcing coefficient (k0 harmonic) numerically from
    # products of order-1 (f) and order-2 (mean h, 2nd-harm g) envelopes.
    # cof(H) entries: cof(H)_11=H22, _12=-H21, _21=-H12, _22=H11.
    # H1 (k0): H11=k0 f1 cos, H12=f1' sin, H21=-k0 f2 sin, H22=f2' cos
    # H2 mean: H11=0, H12=h1', H21=0, H22=h2'   (cos0=1)
    # H2 2nd : H11=2k0 g1 cos2, H12=g1' sin2, H21=-2k0 g2 sin2, H22=g2' cos2
    # p1=f3 cos ; p2_mean=h3 ; p2_2nd=g3 cos2
    # S3 = -(p1 cof(H2) + p2 cof(H1)); we need Div(S3) projected on k0 (cos for comp2 eq, sin for comp1 eq)
    # Assemble Div(S3)_i k0-harmonic via trig product reduction:
    #   cos*1->cos(k0); cos*cos2 -> 1/2[cos(k0)+cos(3k0)]; sin*sin2->1/2[cos(k0)-cos(3k0)]; etc.
    # Keep only k0 terms. Implement S3 components as functions of X2 for the k0 projection.
    # --- S3 tensor, k0-harmonic envelopes (comp,row): we need S3_i1 (mult sin or cos) & S3_i2 ---
    # Compute S3_{iJ} as products, then take Div and project. Use helper for trig reduction.
    # Represent each second-order field's contribution and multiply by order-1.
    # p1*cof(H2): p1=f3 cos(k0)
    #   cof(H2)_11 = H2_22 ; mean h2' (const in X1) and 2nd g2' cos2
    #   product f3 cos * [h2'] -> f3 h2' cos(k0)          (k0)
    #   product f3 cos * [g2' cos2] -> 1/2 f3 g2' cos(k0) (+cos3k0 dropped)
    # so (p1 cofH2)_11 k0-coeff (cos): f3 h2' + 0.5 f3 g2'
    pc11 = f3 * h2p + 0.5 * f3 * g2p
    #   cof(H2)_12 = -H2_21 ; mean 0 ; 2nd -(-2k0 g2 sin2)=2k0 g2 sin2
    #   f3 cos * 2k0 g2 sin2 -> 2k0 f3 g2 * (cos sin2)=1/2[sin3k0 - sin... ] sin-comp k0: cos*sin2=1/2(sin3+sin? ) 
    #   cos a sin2a = 1/2[sin(3a)+sin(a)] -> k0 sin-coeff = 1/2
    pc12 = 0.5 * (2 * k0) * f3 * g2      # (sin k0)
    #   cof(H2)_21 = -H2_12 ; mean -h1' ; 2nd -g1' sin2
    #   f3 cos *(-h1') -> -f3 h1' cos(k0); f3 cos*(-g1' sin2): cos*sin2 sin-coeff 1/2 -> -0.5 f3 g1'
    pc21c = -f3 * h1p                    # (cos k0)
    pc21s = -0.5 * f3 * g1p              # (sin k0)
    #   cof(H2)_22 = H2_11 ; mean 0 ; 2nd 2k0 g1 cos2
    #   f3 cos * 2k0 g1 cos2 -> 0.5*2k0 f3 g1 cos(k0)
    pc22 = 0.5 * (2 * k0) * f3 * g1      # (cos k0)
    # p2*cof(H1): p2_mean=h3 (const), p2_2nd=g3 cos2
    #   cof(H1)_11=H1_22=f2' cos ; h3*f2' cos -> h3 f2' cos(k0); g3 cos2 * f2' cos ->0.5 g3 f2' cos(k0)
    qc11 = h3 * f2p + 0.5 * g3 * f2p
    #   cof(H1)_12=-H1_21=k0 f2 sin ; h3*k0 f2 sin->h3 k0 f2 sin(k0); g3 cos2*k0 f2 sin ->0.5 g3 k0 f2 (sin? cos2 sin=1/2(sin3-sin1)) sin-coeff -0.5
    qc12s = h3 * k0 * f2 - 0.5 * g3 * k0 * f2
    #   cof(H1)_21=-H1_12=-f1' sin ; h3*(-f1' sin)->-h3 f1' sin; g3 cos2*(-f1' sin): cos2 sin sin-coeff -1/2 -> +0.5 g3 f1'
    qc21s = -h3 * f1p + 0.5 * g3 * f1p
    #   cof(H1)_22=H1_11=k0 f1 cos ; h3*k0 f1 cos->h3 k0 f1 cos; g3 cos2*k0 f1 cos->0.5 g3 k0 f1 cos
    qc22 = h3 * k0 * f1 + 0.5 * g3 * k0 * f1
    # S3_{iJ} = -(p1 cofH2 + p2 cofH1)_{iJ}, split into cos(k0)/sin(k0) envelopes
    # row1 (i=1): J=1 uses _11, J=2 uses _12 ; row2 (i=2): _21,_22
    S3_11c = -(pc11 + qc11)              # cos
    S3_12s = -(pc12 + qc12s)             # sin
    S3_21s = -(pc21s + qc21s)            # sin   (pc21 has cos part too -> handle)
    S3_21c = -(pc21c)                    # cos part of S3_21
    S3_22c = -(pc22 + qc22)              # cos
    # Div(S3)_1 = dS3_11/dX1 + dS3_12/dX2 ; comp1 eq carried by sin(k0):
    #   dS3_11(cos)/dX1 = -k0 S3_11c sin ; dS3_12(sin)/dX2 = S3_12s' sin
    F3_1 = -k0 * S3_11c + d(S3_12s)
    # Div(S3)_2 = dS3_21/dX1 + dS3_22/dX2 ; comp2 eq carried by cos(k0):
    #   S3_21 has sin part (->dX1 gives k0 cos) and cos part (->dX1 gives -k0 sin, drop for cos-eq)
    #   dS3_21s(sin)/dX1 = k0 S3_21s cos ; dS3_22(cos)/dX2 = S3_22c' cos
    F3_2 = k0 * S3_21s + d(S3_22c)
    # equilibrium order3: (DivS1(w3))_i = -(DivS3)_i  -> operator RHS = -F3_i
    rI3 = -F3_1
    rII3 = -F3_2
    # constraint order3 cross term: cof(F0):H3 = -(H1:cof? ) -> -(sum of order1xorder2 dets), k0 part
    # det-type cross: contribution to constraint = -[ cof(H1):H2 ]_k0  (since det(F)=1 expands)
    # cof(H1):H2 = H1_22 H2_11 - H1_21 H2_12 - H1_12 H2_21 + H1_11 H2_22 (full contraction of cof)
    # take k0 harmonic of products (mean h and 2nd g):
    cc = (f2p * (0.5 * 2 * k0 * g1)          # H1_22(cos) H2_11(2nd 2k0 g1 cos2): cos*cos2->0.5 cos k0
          + f2p * 0.0                          # mean H2_11=0
          - (-k0 * f2) * (0.5 * g1p)           # -H1_21(sin)*H2_12(2nd g1' sin2): sin*sin2->0.5 cos k0
          - (-k0 * f2) * (h1p)                 # -H1_21(sin? mean H2_12=h1' const): sin*1-> sin, not k0 cos -> drop
          - f1p * (-0.5 * 2 * k0 * g2)         # -H1_12(sin)*H2_21(2nd -2k0 g2 sin2): sin*sin2->0.5cos
          + (k0 * f1) * (0.5 * g2p)            # H1_11(cos)*H2_22(2nd g2' cos2): cos*cos2->0.5cos
          + (k0 * f1) * (h2p))                 # H1_11(cos)*H2_22(mean h2'): cos*1->cos k0
    rIII3 = -cc
    # order3 top BC forcing (gam=0): S1_{i2}(w3) = -(S3_{i2}) at surface
    bcs3 = -S3_12s[0]                          # shear (sin-carried S3_12)
    bcn3 = -S3_22c[0]                          # normal (cos-carried S3_22)

    # --- adjoint & dL/dlam projection ---
    adj = adjoint_mode(N, D, D2, lam, k0, 0.0)
    # numerator: <adj, RHS3> over the assembled RHS vector (interior + BC rows)
    b3 = np.zeros(3 * n)
    b3[0:n] = rI3; b3[n:2 * n] = rII3; b3[2 * n:3 * n] = rIII3
    b3[0] = bcs3; b3[n] = bcn3
    num = adj @ b3
    # denominator: <adj, (dL/dlam) w1>
    e1 = np.concatenate([f1, f2, f3])
    Lp = linear_matrix(N, D, D2, lam + dlam, k0, 0.0)
    Lm = linear_matrix(N, D, D2, lam - dlam, k0, 0.0)
    dLw = ((Lp - Lm) / (2 * dlam)) @ e1
    den = adj @ dLw
    lam2 = -num / den
    return lam2, num, den


if __name__ == "__main__" and False:
    pass
