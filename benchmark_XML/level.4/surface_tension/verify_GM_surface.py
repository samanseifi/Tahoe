"""
Verify Gurtin-Murdoch surface elasticity residual / stiffness against
JAX autodiff for a single 2D edge (2-node line segment).

#54 — Gurtin-Murdoch generalises Young-Laplace by adding a strain-
dependent surface stress.  In 2D the surface is a 1D curve, so the
surface strain is a scalar and one effective surface stiffness suffices:

    E_s := 2 mu_s + lambda_s        (units: N/m)

Reference and current edge lengths:
    L_0 = || X_2 - X_1 ||
    L   = || (X_2 + u_2) - (X_1 + u_1) ||

Engineering surface strain:
    eps_eng = (L - L_0) / L_0

Surface free energy (per unit reference length, integrated over L_0):
    W_s = gamma_0 * L  +  (1/2) * E_s * (L - L_0)^2 / L_0
        = L_0 * [ gamma_0 * (1 + eps_eng)  +  (1/2) * E_s * eps_eng^2 ]

The first term is the classical Young-Laplace contribution (recovers the
existing implementation when E_s = 0).  The second term is the new GM
elastic stretch energy.

Surface (1st PK) stress in the surface (force per unit reference length):
    sigma_s = dW_s/dL = gamma_0 + E_s * eps_eng

Nodal residual contributed by the surface energy (4 DOFs INTERLEAVED:
[u_x1, u_y1, u_x2, u_y2]):
    R_i = dW_s/du_i = sigma_s * dL/du_i

with dL/du_i:
    dL/du_x1 = (x1 - x2) / L
    dL/du_y1 = (y1 - y2) / L
    dL/du_x2 = (x2 - x1) / L
    dL/du_y2 = (y2 - y1) / L

Surface tangent stiffness:
    K_ij = d^2 W_s / du_i du_j
         = (E_s / L_0) * (dL/du_i)(dL/du_j)       <-- NEW: material term
           + sigma_s * d^2 L / du_i du_j           <-- geometric, sigma_s replaces gamma

In terms of the existing code's fB = L * dL/du and the K1 = identity-like
matrix already encoded:
    K_geom     = (sigma_s / L) * K1 + (-sigma_s / L^3) * fB outer fB
    K_material = (E_s / (L_0 * L^2)) * fB outer fB
    K_total    = K_geom + K_material

In the Young-Laplace limit (E_s = 0, sigma_s = gamma_0), this collapses
to the existing manual_K() from verify_surface_tension.py exactly.

DOF sign convention: SimoQ1P0_Surface stores fD = -dW_s/du; fRHS += fD.
This script verifies both the residual and the stiffness against jax
autodiff so we can lock in the GM extension before writing C++.
"""

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)


# ── 1. Energy ────────────────────────────────────────────────────────────────
def Wsurface(u, X, gamma_0, E_s):
    """Gurtin-Murdoch surface energy on a 2-node edge."""
    x1 = X[0] + u[0]; y1 = X[1] + u[1]
    x2 = X[2] + u[2]; y2 = X[3] + u[3]
    L   = jnp.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    L_0 = jnp.sqrt((X[0] - X[2]) ** 2 + (X[1] - X[3]) ** 2)
    eps_eng = (L - L_0) / L_0
    return gamma_0 * L + 0.5 * E_s * (L - L_0) ** 2 / L_0


# ── 2. Closed-form residual and stiffness (what the C++ should code) ────────
def manual_fD(u, X, gamma_0, E_s):
    """Residual fD = -dW_s/du (added to fRHS, matches existing convention)."""
    x1 = X[0] + u[0]; y1 = X[1] + u[1]
    x2 = X[2] + u[2]; y2 = X[3] + u[3]
    L   = jnp.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    L_0 = jnp.sqrt((X[0] - X[2]) ** 2 + (X[1] - X[3]) ** 2)
    eps_eng = (L - L_0) / L_0
    sigma_s = gamma_0 + E_s * eps_eng            # NEW: strain-dependent surface stress
    c = -sigma_s / L                              # was -gamma_0/L in Young-Laplace
    return jnp.array([
        c * (x1 - x2),
        c * (y1 - y2),
        c * (x2 - x1),
        c * (y2 - y1),
    ])


def manual_K(u, X, gamma_0, E_s):
    """Stiffness K_total = K_geom + K_material."""
    x1 = X[0] + u[0]; y1 = X[1] + u[1]
    x2 = X[2] + u[2]; y2 = X[3] + u[3]
    L   = jnp.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    L_0 = jnp.sqrt((X[0] - X[2]) ** 2 + (X[1] - X[3]) ** 2)
    eps_eng = (L - L_0) / L_0
    sigma_s = gamma_0 + E_s * eps_eng

    K1 = jnp.array([
        [ 1.,  0., -1.,  0.],
        [ 0.,  1.,  0., -1.],
        [-1.,  0.,  1.,  0.],
        [ 0., -1.,  0.,  1.],
    ])
    fB     = jnp.array([x1 - x2, y1 - y2, x2 - x1, y2 - y1])
    fB_outer = jnp.outer(fB, fB)

    K_geom     = (sigma_s / L) * K1 - (sigma_s / L ** 3) * fB_outer
    K_material = (E_s / (L_0 * L ** 2)) * fB_outer
    return K_geom + K_material


# ── 3. Run tests ─────────────────────────────────────────────────────────────
def run_one(label, X, u, gamma_0, E_s, atol=1e-12):
    grad_jax = jax.grad(Wsurface)(u, X, gamma_0, E_s)
    Hess_jax = jax.hessian(Wsurface)(u, X, gamma_0, E_s)

    fD  = manual_fD(u, X, gamma_0, E_s)
    K   = manual_K (u, X, gamma_0, E_s)

    err_res = jnp.max(jnp.abs(grad_jax + fD))        # fD == -grad
    err_K   = jnp.max(jnp.abs(Hess_jax - K))

    sym_err = jnp.max(jnp.abs(K - K.T))               # tangent symmetry

    pass_res = bool(err_res < atol)
    pass_K   = bool(err_K   < atol)
    pass_sym = bool(sym_err < 1e-14)

    tag = "PASS" if (pass_res and pass_K and pass_sym) else "FAIL"
    print(f"  [{tag}] {label:55s} "
          f"res err = {err_res:.2e}  K err = {err_K:.2e}  sym err = {sym_err:.2e}")

    if not (pass_res and pass_K and pass_sym):
        print("    JAX grad:   ", np.array(grad_jax))
        print("    manual fD:  ", np.array(fD))
        print("    JAX Hess:\n", np.array(Hess_jax))
        print("    manual K:\n", np.array(K))

    return pass_res and pass_K and pass_sym


def run_all():
    print("=" * 78)
    print("Gurtin-Murdoch surface verification (#54) — JAX autodiff vs manual")
    print("DOF ordering: INTERLEAVED  [u_x1, u_y1, u_x2, u_y2]")
    print("=" * 78)

    cases = []
    # ── Reference (Young-Laplace) cases: E_s = 0 must match existing impl ──
    for label, X, u in [
        ("YL  horizontal edge, rest",        [0., 0., 1., 0.], [0., 0., 0., 0.]),
        ("YL  horizontal edge, perturbed",   [0., 0., 1., 0.], [0., 0.1, 0., -0.05]),
        ("YL  diagonal edge",                [0., 0., 0.6, 0.8], [0.05, -0.03, -0.02, 0.04]),
        ("YL  large deformation",            [0., 0., 1., 0.], [-0.3, 0.4, 0.2, -0.1]),
    ]:
        cases.append((label, jnp.array(X), jnp.array(u), 0.05, 0.0))

    # ── GM-specific cases: E_s > 0 stiff surface ────────────────────────────
    for label, X, u, gamma_0, E_s in [
        ("GM  horizontal rest, mu_s+lam_s=1",  [0., 0., 1., 0.], [0., 0., 0., 0.], 0.05, 1.0),
        ("GM  horizontal +5% stretch, E_s=1",  [0., 0., 1., 0.], [-0.025, 0., 0.025, 0.],
         0.05, 1.0),
        ("GM  horizontal -5% compress, E_s=2", [0., 0., 1., 0.], [0.025, 0., -0.025, 0.],
         0.05, 2.0),
        ("GM  diagonal +10%, gamma=2 E_s=5",   [0., 0., 0.6, 0.8], [-0.03, -0.04, 0.03, 0.04],
         2.0, 5.0),
        ("GM  large twist, E_s=10",            [0., 0., 1., 0.], [-0.2, 0.4, 0.2, -0.3],
         0.5, 10.0),
        ("GM  negative E_s (softening)",       [0., 0., 1., 0.], [-0.05, 0.0, 0.05, 0.0],
         0.05, -0.5),
        ("GM  gamma=0 pure surface elasticity",[0., 0., 1., 0.], [-0.025, 0., 0.025, 0.],
         0.0, 1.0),
    ]:
        cases.append((label, jnp.array(X), jnp.array(u), gamma_0, E_s))

    print("\n--- Suite ---")
    n_pass = 0
    for label, X, u, gamma_0, E_s in cases:
        if run_one(label, X, u, gamma_0, E_s):
            n_pass += 1

    print()
    print(f"Result: {n_pass}/{len(cases)} passed")
    print("=" * 78)

    # ── Sanity prints ──────────────────────────────────────────────────────
    print("\nManual K at horizontal rest, gamma_0=0.05, E_s=1.0:")
    X0 = jnp.array([0., 0., 1., 0.]); u0 = jnp.zeros(4)
    K0 = manual_K(u0, X0, 0.05, 1.0)
    print(np.array(K0))
    print("  Diagonal:", np.diag(np.array(K0)))
    print("  Expected non-zero diagonal terms include:")
    print("    K_xx  =  E_s / L_0    = 1.0      (material stretch stiffness in x)")
    print("    K_yy  =  gamma_0 / L  = 0.05     (Young-Laplace transverse restoring)")

    print("\nGM limit check (E_s -> 0):")
    K_YL = manual_K(u0, X0, 0.05, 0.0)
    print("  K with E_s=0:")
    print(np.array(K_YL))
    print("  -> must equal the existing verify_surface_tension.py manual_K() bitwise.")

    print("\nSurface stress samples sigma_s = gamma_0 + E_s * eps_eng:")
    for stretch in [0.9, 1.0, 1.1, 1.2]:
        eps = stretch - 1.0
        sig = 0.05 + 1.0 * eps
        print(f"  L/L_0={stretch:.2f}, eps={eps:+.3f}  ->  sigma_s = {sig:+.3f}")
    print("  At L/L_0=1: sigma_s = gamma_0 (Young-Laplace recovered).")
    print("  At L/L_0=1.1, E_s=1: sigma_s = 0.15 (3x the bare gamma_0).")


if __name__ == "__main__":
    run_all()
