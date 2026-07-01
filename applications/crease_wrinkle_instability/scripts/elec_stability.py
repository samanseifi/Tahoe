#!/usr/bin/env python3
r"""
Electromechanical surface-WRINKLING of a pre-strained incompressible neo-Hookean
dielectric-elastomer film with surface tension, for finite in-plane pre-stretch
lambda.

Model (reference configuration, film 0<=X2<=H, mu=eps=1):
  * the field is divergence-free in the bulk, so the perturbation obeys the
    purely mechanical envelope equation with roots q in {+-K, +-K/lam^2};
    f2(X2)=sum c_i e^{q_i X2}, f1=-(lam^2/K) f2'.
  * mechanical nominal tractions on a horizontal face:
        t1(shear)  = -[(lam^2/K) f2'' + (K/lam^2) f2]
        t2(normal) =  (2+lam^4) f2' - (lam^4/K^2) f2'''
  * boundary conditions
        base  X2=H (bonded/clamped):  f2 = 0,  f2' = 0
        top   X2=0 (free, electroded): t1 = 0  (shear-free)
              t2 + (eps E2^2 lam^2/H) f2 - (gamma K^2/lam) f2 = 0
    where E2 = true (current) field = lam*Etilde, Etilde=Phi/H the nominal field,
    the Maxwell term  eps E2^2 lam^2/H  is the gap-change electrostatic traction,
    and gamma K^2/lam is the (stabilizing) Young-Laplace traction, gamma=gbar*mu*H.

The field enters the top normal condition linearly through s = eps E2^2/mu, so for
each K the critical s solves a 4x4 determinant linear in s; the wrinkle threshold
is the minimum over K.  At lam=1, gbar=0/0.5/1/2/5/10/20 it gives
  Ec*sqrt(eps/mu) = 2.50, 2.79, 2.96, 3.19, 3.62, 4.04, 4.57   (grows as sqrt(gbar)).

Conventions reported:  true field E2c*sqrt(eps/mu)=sqrt(s_c);
nominal field  Etilde_c*sqrt(eps/mu) = E2c*sqrt(eps/mu)/lam.
"""
import argparse
import numpy as np


def _Mat(K, lam, gbar, H, s):
    """4x4 boundary matrix (mu=1); columns are modes q in {K,-K,K/lam^2,-K/lam^2}."""
    qs = (K, -K, K / lam**2, -K / lam**2)
    gamma = gbar * H
    M = np.zeros((4, 4))
    for j, q in enumerate(qs):
        eqH = np.exp(min(q * H, 700.0))
        M[0, j] = eqH                                   # base f2(H)=0
        M[1, j] = q * eqH                               # base f2'(H)=0
        M[2, j] = -((lam**2 / K) * q**2 + (K / lam**2))  # top shear-free t1(0)=0
        t2 = (2 + lam**4) * q - (lam**4 / K**2) * q**3
        M[3, j] = t2 + (s * lam**2 / H) - (gamma * K**2 / lam)  # top normal
    return M


def _s_crit_at_K(K, lam, gbar, H):
    """critical s = eps E2^2/mu at fixed K (det is linear in s)."""
    M0 = _Mat(K, lam, gbar, H, 0.0)
    M1 = _Mat(K, lam, gbar, H, 1.0)
    sc = np.max(np.abs(M0), axis=0)
    sc[sc == 0] = 1.0
    with np.errstate(over="ignore", invalid="ignore"):
        d0 = np.linalg.det(M0 / sc)
        d1 = np.linalg.det(M1 / sc)
    coeff = d1 - d0
    if not np.isfinite(coeff) or abs(coeff) < 1e-300:
        return None
    s = -d0 / coeff
    return s if (np.isfinite(s) and s > 0) else None


def _nudge(lam):
    return 1.0 + 1e-4 if abs(lam - 1.0) < 1e-4 else lam


def true_field_at_K(lam, K, gbar, H):
    """neutral true field E2*sqrt(eps/mu)=sqrt(s_c) at wavenumber K (or None)."""
    s = _s_crit_at_K(K, _nudge(lam), gbar, H)
    return np.sqrt(s) if s is not None else None


def neutral_curve(lam, gbar, H, Ks):
    """(wavelength l/H_f, nominal field Etilde*sqrt(eps/mu)) along the neutral curve.
    l/H_f is the CURRENT (deformed) wavelength: the in-plane stretch is fixed at the
    pre-strain lam=lambda_pre (the voltage only thins the film), so the deformed
    crest-to-crest spacing is lam*(2*pi/K), normalized by the reference H_f."""
    lam = _nudge(lam)
    ll, EE = [], []
    for K in Ks:
        E2 = true_field_at_K(lam, K, gbar, H)
        if E2 is not None:
            ll.append(2 * np.pi / K / H * lam)   # current (deformed) wavelength / H_f
            EE.append(E2 / lam)                  # nominal = true / lam
    return np.array(ll), np.array(EE)


def critical_field(lam, gbar, H, Ks):
    """critical NOMINAL field Etilde_c*sqrt(eps/mu) (min over K) and selected K_c."""
    lam = _nudge(lam)
    best = None
    for K in Ks:
        s = _s_crit_at_K(K, lam, gbar, H)
        if s is not None and (best is None or s < best[0]):
            best = (s, K)
    if best is None:
        return None
    return np.sqrt(best[0]) / lam, best[1]       # nominal field, K_c


def critical_true_field(lam, gbar, H, Ks):
    """critical TRUE field E2c*sqrt(eps/mu)=sqrt(s_c) (min over K) and K_c."""
    r = critical_field(lam, gbar, H, Ks)
    return (r[0] * _nudge(lam), r[1]) if r else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--gbar", type=float, default=0.0)
    ap.add_argument("--H", type=float, default=4.0)
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()
    Ks = np.linspace(0.02, 3.0, 1500)

    if args.validate:
        print("== lambda=1 wrinkle threshold Ec*sqrt(eps/mu) vs gbar ==")
        print("   gbar :    Ec    l/H")
        for gb in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
            r = critical_true_field(1.0, gb, args.H, Ks)
            print(f"   {gb:>5}: {r[0]:>6.3f}  {2*np.pi/r[1]/args.H:>5.2f}")
        print("== pre-strain trend at gbar=1 (true E2c and nominal Etilde_c) ==")
        for lam in (0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20):
            r = critical_field(lam, 1.0, args.H, Ks)
            if r:
                print(f"   lam={lam:.2f} (eps_pre={1-lam:+.2f}): E2c={r[0]*_nudge(lam):.3f}  "
                      f"Etilde_c={r[0]:.3f}  l/H={2*np.pi/r[1]/args.H:.2f}")
        return

    r = critical_field(args.lam, args.gbar, args.H, Ks)
    if r:
        E2 = r[0] * _nudge(args.lam)
        print(f"lam={args.lam} gbar={args.gbar}: true E2c*sqrt(eps/mu)={E2:.3f}, "
              f"nominal Etilde_c*sqrt(eps/mu)={r[0]:.3f}, K_c={r[1]:.2f}, "
              f"l/H={2*np.pi/r[1]/args.H:.2f}")
    else:
        print("no bifurcation found")


if __name__ == "__main__":
    main()
