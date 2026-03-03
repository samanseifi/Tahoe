"""
Verify surface tension residual (FormKd) and stiffness (FormStiffness)
using JAX automatic differentiation.

Energy: Psi = gamma * L
where L = sqrt((x1-x2)^2 + (y1-y2)^2) is the current edge length,
and xi = Xi + ui (current = reference + displacement).

DOF vector uses INTERLEAVED ordering (confirmed from FieldT::SetLocalEqnos):
  u = [u_x1, u_y1, u_x2, u_y2]   (4 DOFs for 2-node 2D edge)

The residual contribution R = dPsi/du  (internal force, positive = pulls inward)
The stiffness           K = d^2Psi/du^2

We compare these against the manual formulas coded in SimoQ1P0_Surface:
  fD[0] = -gamma/L * (x1 - x2)   (code uses coeff3 = -gamma/L)
  fD[1] = -gamma/L * (y1 - y2)
  fD[2] = -gamma/L * (x2 - x1)
  fD[3] = -gamma/L * (y2 - y1)

  K_code = (gamma/L) * K1 + (-gamma/L^3) * (fB outer fB)
  where K1 = [[1,0,-1,0],[0,1,0,-1],[-1,0,1,0],[0,-1,0,1]]
  and   fB = [x1-x2, y1-y2, x2-x1, y2-y1]

NOTE: fRHS accumulates -dPsi/du (the force as seen from the solver perspective),
so fD = -dPsi/du. We verify that below.
"""

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)


# ── energy function ──────────────────────────────────────────────────────────
def psi(u, X, gamma):
    """
    Surface energy Psi = gamma * L_current.

    Parameters
    ----------
    u : (4,) array  — displacements [u_x1, u_y1, u_x2, u_y2] (INTERLEAVED)
    X : (4,) array  — reference coords [X1, Y1, X2, Y2]
    gamma : float
    """
    x1 = X[0] + u[0]
    y1 = X[1] + u[1]
    x2 = X[2] + u[2]
    y2 = X[3] + u[3]
    L  = jnp.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    return gamma * L


# ── manual formulas from SimoQ1P0_Surface ───────────────────────────────────
def manual_fD(u, X, gamma):
    """
    Manual residual vector fD as coded (= -dPsi/du).
    fRHS += fD  →  fRHS gets -dPsi/du  (internal force sign convention).
    """
    x1 = X[0] + u[0]; y1 = X[1] + u[1]
    x2 = X[2] + u[2]; y2 = X[3] + u[3]
    L  = jnp.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    c  = -gamma / L
    return jnp.array([c*(x1-x2), c*(y1-y2), c*(x2-x1), c*(y2-y1)])


def manual_K(u, X, gamma):
    """
    Manual stiffness K_Total as coded.
    K_Total = (gamma/L)*K1 + (-gamma/L^3)*(fB outer fB)
    K1 = [[1,0,-1,0],[0,1,0,-1],[-1,0,1,0],[0,-1,0,1]]
    """
    x1 = X[0] + u[0]; y1 = X[1] + u[1]
    x2 = X[2] + u[2]; y2 = X[3] + u[3]
    L  = jnp.sqrt((x1-x2)**2 + (y1-y2)**2)

    K1 = jnp.array([
        [ 1.,  0., -1.,  0.],
        [ 0.,  1.,  0., -1.],
        [-1.,  0.,  1.,  0.],
        [ 0., -1.,  0.,  1.],
    ])
    fB = jnp.array([x1-x2, y1-y2, x2-x1, y2-y1])
    K2 = jnp.outer(fB, fB)

    return (gamma / L) * K1 + (-gamma / L**3) * K2


# ── run tests ────────────────────────────────────────────────────────────────
def run_tests():
    gamma = 0.05

    test_cases = [
        # (label, X_ref, u_disp)
        ("horizontal edge, no disp",
         jnp.array([0., 0., 1., 0.]),
         jnp.array([0., 0., 0., 0.])),

        ("horizontal edge, y-perturbed",
         jnp.array([0., 0., 1., 0.]),
         jnp.array([0., 0.1, 0., -0.05])),

        ("diagonal edge",
         jnp.array([0., 0., 0.6, 0.8]),
         jnp.array([0.05, -0.03, -0.02, 0.04])),

        ("large deformation",
         jnp.array([0., 0., 1., 0.]),
         jnp.array([-0.3, 0.4, 0.2, -0.1])),
    ]

    print("=" * 70)
    print("Surface tension verification: JAX autodiff vs manual formulas")
    print("DOF ordering: INTERLEAVED  [u_x1, u_y1, u_x2, u_y2]")
    print("=" * 70)

    all_passed = True

    for label, X, u in test_cases:
        print(f"\n--- {label} ---")

        # JAX gradient  =  dPsi/du  (this is the +internal force direction)
        jax_grad = jax.grad(psi)(u, X, gamma)
        # Manual fD = -dPsi/du  (what the code adds to fRHS)
        man_fD   = manual_fD(u, X, gamma)

        grad_err = jnp.max(jnp.abs(jax_grad + man_fD))   # fD == -grad
        print(f"  Residual check  (fD == -dPsi/du):  max|fD + grad| = {grad_err:.3e}", end="")

        # JAX Hessian  =  d^2Psi/du^2
        jax_H    = jax.hessian(psi)(u, X, gamma)
        # Manual K_Total  =  d^2Psi/du^2
        man_K    = manual_K(u, X, gamma)

        K_err = jnp.max(jnp.abs(jax_H - man_K))
        print(f"   |  Stiffness check (K == d²Psi/du²): max|K - H| = {K_err:.3e}", end="")

        ok_r = float(grad_err) < 1e-12
        ok_k = float(K_err)    < 1e-12
        status = "PASS" if (ok_r and ok_k) else "FAIL"
        print(f"   [{status}]")

        if not (ok_r and ok_k):
            all_passed = False
            print(f"\n  JAX grad:  {np.array(jax_grad)}")
            print(f"  man fD:    {np.array(man_fD)}")
            print(f"\n  JAX Hess:\n{np.array(jax_H)}")
            print(f"  man K:\n{np.array(man_K)}")

    print("\n" + "=" * 70)
    print("Overall:", "ALL PASS" if all_passed else "SOME FAILURES")
    print("=" * 70)

    # ── extra: show what K looks like for a horizontal edge at rest ──────
    print("\nK_Total for horizontal edge (X=[0,0,1,0], u=0):")
    X0 = jnp.array([0., 0., 1., 0.])
    u0 = jnp.zeros(4)
    K0 = manual_K(u0, X0, gamma)
    H0 = jax.hessian(psi)(u0, X0, gamma)
    print(f"  Manual K:\n{np.array(K0)}")
    print(f"  JAX H  :\n{np.array(H0)}")
    print(f"  Diagonal of K: {np.diag(np.array(K0))}")
    print("  → In-plane (x) stiffness is zero (correct: moving along edge")
    print("    doesn't change length to first order when edge is straight).")
    print("  → Transverse (y) stiffness is gamma/L = curvature restoring.")

    # ── consistency check: K should equal d(fRHS)/du sign-wise ─────────
    print("\nConsistency: K = +d(dPsi/du)/du  → K*du = d(fRHS)/du * (-1) * (-du) = K*du")
    print("  (FormStiffness adds K positively to fLHS. Solver uses K*Δu = -fRHS.)")
    print("  d(fRHS)/du = d(-dPsi/du)/du = -d²Psi/du² = -K")
    print("  So adding K to fLHS is CORRECT: residual decreases when K*Δu = -fRHS.")


if __name__ == "__main__":
    run_tests()
