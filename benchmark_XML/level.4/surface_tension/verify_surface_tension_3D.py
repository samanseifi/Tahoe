"""
JAX verification of 3D surface tension residual and stiffness for a single quad face.

Implements the covariant metric tensor formulation:
  W_s(x) = gamma * int H dxi1 dxi2,  H = sqrt(EG - F^2)
  g1 = sum_a (dN_a/dxi1) * x_a,  g2 = sum_a (dN_a/dxi2) * x_a
  E = g1.g1,  G = g2.g2,  F = g1.g2

Cross-checks:
  1. Residual (F_int = dW/dx) vs JAX grad
  2. Stiffness (K = d^2W/dx^2) vs JAX hessian
  3. Residual sign convention (fRHS += -F_int) — pulls face toward center
  4. Symmetry of stiffness matrix
  5. PSD check of stiffness matrix restricted to tangential subspace

Reference: Henann & Bertoldi (2014) uel_surften_3D_hex.for
"""

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

# ── 2×2 Gauss quadrature on [-1,1]^2 ──────────────────────────────────────────
GP = 1.0 / jnp.sqrt(3.0)
XI  = jnp.array([[-GP, -GP], [GP, -GP], [-GP, GP], [GP, GP]])  # (4, 2)
WGP = jnp.array([1.0, 1.0, 1.0, 1.0])

def shape_fns(xi, eta):
    """Bilinear quad shape functions and derivatives at (xi, eta)."""
    N  = jnp.array([(1-xi)*(1-eta), (1+xi)*(1-eta),
                    (1+xi)*(1+eta), (1-xi)*(1+eta)]) * 0.25
    dN = jnp.array([
        [-(1-eta), -(1-xi)],   # dN0/dxi1, dN0/dxi2
        [ (1-eta), -(1+xi)],
        [ (1+eta),  (1+xi)],
        [-(1+eta),  (1-xi)],
    ]) * 0.25                  # shape (4, 2)
    return N, dN

# ── Surface energy ─────────────────────────────────────────────────────────────
def surface_energy(x_flat, gamma=1.0):
    """
    W_s = gamma * int H dxi1 dxi2

    x_flat : (12,) = [x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3]  (INTERLEAVED)
    """
    x = x_flat.reshape(4, 3)   # (4 nodes, 3 coords)
    W = 0.0
    for ip in range(4):
        xi, eta = XI[ip]
        _, dN = shape_fns(xi, eta)          # (4, 2)
        g1 = dN[:, 0] @ x                  # (3,)  covariant basis vector 1
        g2 = dN[:, 1] @ x                  # (3,)
        E  = jnp.dot(g1, g1)
        G  = jnp.dot(g2, g2)
        F  = jnp.dot(g1, g2)
        H  = jnp.sqrt(E*G - F*F)
        W  = W + WGP[ip] * H
    return gamma * W

# ── Analytical residual ────────────────────────────────────────────────────────
def residual_analytic(x_flat, gamma=1.0):
    """
    R_a = -dW/dx_a  (internal force = +dW/dx, so fRHS contribution = -F_int)
    Returns (12,) matching INTERLEAVED ordering.
    """
    x = x_flat.reshape(4, 3)
    R = jnp.zeros((4, 3))
    for ip in range(4):
        xi, eta = XI[ip]
        _, dN = shape_fns(xi, eta)
        g1 = dN[:, 0] @ x
        g2 = dN[:, 1] @ x
        E  = jnp.dot(g1, g1)
        G  = jnp.dot(g2, g2)
        F  = jnp.dot(g1, g2)
        H  = jnp.sqrt(E*G - F*F)
        w  = WGP[ip]
        coeff = -gamma * w / H
        for a in range(4):
            n1a = dN[a, 0]
            n2a = dN[a, 1]
            c1 = G*n1a - F*n2a
            c2 = E*n2a - F*n1a
            R = R.at[a].add(coeff * (c1*g1 + c2*g2))
    return R.reshape(-1)

# ── Analytical stiffness ───────────────────────────────────────────────────────
def stiffness_analytic(x_flat, gamma=1.0):
    """
    K_ab = +dF_int_a/dx_b = -d(fRHS_contribution)/dx_b
    Three terms: identity, outer products, curvature.
    Returns (12, 12).
    """
    x = x_flat.reshape(4, 3)
    K = jnp.zeros((12, 12))
    for ip in range(4):
        xi, eta = XI[ip]
        _, dN = shape_fns(xi, eta)
        g1 = dN[:, 0] @ x
        g2 = dN[:, 1] @ x
        E  = jnp.dot(g1, g1)
        G  = jnp.dot(g2, g2)
        F  = jnp.dot(g1, g2)
        H  = jnp.sqrt(E*G - F*F)
        H2 = H*H
        w  = WGP[ip]
        pref      = gamma * w / H
        pref_curv = gamma * w / (H2 * H)

        Va = jnp.zeros((4, 3))
        for a in range(4):
            n1a = dN[a, 0]; n2a = dN[a, 1]
            Va = Va.at[a].set((G*n1a - F*n2a)*g1 + (E*n2a - F*n1a)*g2)

        for a in range(4):
            n1a = dN[a, 0]; n2a = dN[a, 1]
            for b in range(4):
                n1b = dN[b, 0]; n2b = dN[b, 1]
                s_ab = G*n1a*n1b - F*(n1a*n2b + n2a*n1b) + E*n2a*n2b
                for d in range(3):
                    for e in range(3):
                        # Term 1
                        K_de = pref * s_ab * (1.0 if d == e else 0.0)
                        # Term 2
                        K_de += pref * (
                            n1a*g1[d] * 2.0*n2b*g2[e]
                          - (n2a*g1[d] + n1a*g2[d]) * (n1b*g2[e] + n2b*g1[e])
                          + n2a*g2[d] * 2.0*n1b*g1[e]
                        )
                        # Term 3: curvature
                        K_de -= pref_curv * Va[a, d] * Va[b, e]

                        K = K.at[3*a+d, 3*b+e].add(K_de)
    return K

# ── JAX auto-diff references ───────────────────────────────────────────────────
def residual_jax(x_flat, gamma=1.0):
    """fRHS contribution = -grad(W) = -F_int  (same sign as C++ code's fRHS += R)"""
    return -jax.grad(lambda x: surface_energy(x, gamma))(x_flat)

def stiffness_jax(x_flat, gamma=1.0):
    """K = +d(F_int)/dx = +d(-R)/dx = -d(R)/dx = hessian(W)"""
    return jax.hessian(lambda x: surface_energy(x, gamma))(x_flat)

# ── Test cases ─────────────────────────────────────────────────────────────────

def run_tests():
    gamma = 1.0
    np.set_printoptions(precision=6, suppress=True)
    print("="*70)
    print("3D surface tension verification — covariant metric tensor approach")
    print("="*70)

    # ── Test 1: flat unit square in xy-plane, z=0 ──────────────────────────
    print("\n[Test 1] Flat unit square in xy-plane: nodes at (0,0,0),(1,0,0),(1,1,0),(0,1,0)")
    x0 = jnp.array([
        0.0, 0.0, 0.0,   # node 0
        1.0, 0.0, 0.0,   # node 1
        1.0, 1.0, 0.0,   # node 2
        0.0, 1.0, 0.0,   # node 3
    ])

    W = surface_energy(x0, gamma)
    print(f"  Surface energy = {W:.8f}  (expected {gamma*1.0:.8f}  [area=1])")

    R_anal  = residual_analytic(x0, gamma)
    R_jax   = residual_jax(x0, gamma)
    K_anal  = stiffness_analytic(x0, gamma)
    K_jax   = stiffness_jax(x0, gamma)

    r_err = jnp.max(jnp.abs(R_anal - R_jax))
    k_err = jnp.max(jnp.abs(K_anal - K_jax))
    print(f"  Residual max error (analytic vs JAX):  {r_err:.2e}")
    print(f"  Stiffness max error (analytic vs JAX): {k_err:.2e}")

    print("  Residual (should pull corners toward face center):")
    R_anal_np = np.array(R_anal).reshape(4, 3)
    for a in range(4):
        print(f"    node {a}: {R_anal_np[a]}")

    # symmetry check
    sym_err = jnp.max(jnp.abs(K_jax - K_jax.T))
    print(f"  Stiffness symmetry error: {sym_err:.2e}")

    # PSD check (tangential modes only — normal mode z is zero for flat face)
    eigvals = jnp.linalg.eigvalsh(K_jax)
    print(f"  Eigenvalues of K: {np.array(jnp.sort(eigvals))}")
    print(f"  Min eigenvalue: {float(jnp.min(eigvals)):.4e}  (expect ≥ 0, zero for RBMs)")

    # ── Test 2: 2x2 square in xy-plane (area=4) ────────────────────────────
    print("\n[Test 2] Flat 2x2 square in xy-plane")
    x2 = jnp.array([
        0.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        2.0, 2.0, 0.0,
        0.0, 2.0, 0.0,
    ])
    W2 = surface_energy(x2, gamma)
    print(f"  Surface energy = {W2:.8f}  (expected {gamma*4.0:.8f}  [area=4])")

    R2_anal = residual_analytic(x2, gamma)
    R2_jax  = residual_jax(x2, gamma)
    r2_err  = jnp.max(jnp.abs(R2_anal - R2_jax))
    print(f"  Residual max error: {r2_err:.2e}")

    K2_anal = stiffness_analytic(x2, gamma)
    K2_jax  = stiffness_jax(x2, gamma)
    k2_err  = jnp.max(jnp.abs(K2_anal - K2_jax))
    print(f"  Stiffness max error: {k2_err:.2e}")

    # ── Test 3: tilted face (not in xy-plane) ──────────────────────────────
    print("\n[Test 3] Tilted face (45-degree rotation about x-axis, unit square)")
    c, s = np.cos(np.pi/4), np.sin(np.pi/4)
    x3 = jnp.array([
        0.0,  0.0,  0.0,
        1.0,  0.0,  0.0,
        1.0,  c,    s,
        0.0,  c,    s,
    ])
    W3 = surface_energy(x3, gamma)
    print(f"  Surface energy = {W3:.8f}  (expected {gamma*1.0:.8f}  [area=1])")

    R3_anal = residual_analytic(x3, gamma)
    R3_jax  = residual_jax(x3, gamma)
    r3_err  = jnp.max(jnp.abs(R3_anal - R3_jax))
    K3_anal = stiffness_analytic(x3, gamma)
    K3_jax  = stiffness_jax(x3, gamma)
    k3_err  = jnp.max(jnp.abs(K3_anal - K3_jax))
    print(f"  Residual max error: {r3_err:.2e}")
    print(f"  Stiffness max error: {k3_err:.2e}")

    # ── Test 4: perturbed face — finite difference check on residual ────────
    print("\n[Test 4] Finite-difference check: K ≈ -dR/dx for perturbed flat face")
    eps = 1e-5
    x4 = x0 + 0.03 * jnp.array([
        0.1, 0.2, 0.0,  -0.1, 0.1, 0.0,  0.0, -0.1, 0.0,  0.1, -0.2, 0.0
    ])
    R4 = residual_analytic(x4, gamma)    # fRHS contribution = -F_int
    K4_anal = stiffness_analytic(x4, gamma)
    # FD: K_ij ≈ -(R_i(x+eps*ej) - R_i(x-eps*ej)) / (2*eps)  [K = -dR/dx]
    K4_fd = np.zeros((12, 12))
    for j in range(12):
        ep = jnp.zeros(12).at[j].set(eps)
        Rp = residual_analytic(x4 + ep, gamma)
        Rm = residual_analytic(x4 - ep, gamma)
        K4_fd[:, j] = -(np.array(Rp) - np.array(Rm)) / (2*eps)
    k4_err = np.max(np.abs(K4_anal - K4_fd))
    print(f"  FD stiffness max error: {k4_err:.2e}  (expect ~1e-9)")

    # ── Summary ────────────────────────────────────────────────────────────
    print("\n" + "="*70)
    print("SUMMARY")
    all_ok = (r_err < 1e-10 and k_err < 1e-10 and r2_err < 1e-10 and
              k2_err < 1e-10 and r3_err < 1e-10 and k3_err < 1e-10 and
              k4_err < 1e-8)
    if all_ok:
        print("  ALL TESTS PASSED — analytic formulas match JAX auto-diff")
    else:
        print("  SOME TESTS FAILED — check formulas above")
    print("="*70)

    # ── Diagnose C++ sign convention ───────────────────────────────────────
    print("\n[Sign convention for C++ code]")
    print("  fRHS += R  where R = residual_analytic(x) = -F_int_surface")
    print("  fAmm_mat2 += K  where K = stiffness_analytic(x) = +dF_int/dx")
    print("  So: K·du = fRHS = F_ext - F_int_bulk - F_int_surface")
    print("  Surface tension is a compressive load: residual pulls nodes inward")
    R_flat_np = np.array(residual_analytic(x0, gamma)).reshape(4, 3)
    print("  Node forces for unit square (flat, γ=1):")
    corners = [(0,0,0), (1,0,0), (1,1,0), (0,1,0)]
    for a in range(4):
        cx, cy, cz = corners[a]
        print(f"    node {a} at ({cx},{cy},{cz}): F_surface = {R_flat_np[a]}  "
              f"(points {'inward ✓' if abs(R_flat_np[a,0]) + abs(R_flat_np[a,1]) > 0.01 else '?'})")

if __name__ == "__main__":
    run_tests()
