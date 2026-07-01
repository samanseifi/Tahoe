#!/usr/bin/env python3
r"""
Direct nonlinear continuation of the surface-wrinkling branch of a compressed
incompressible neo-Hookean film with surface tension.  Robust alternative to the
asymptotic Landau coefficient: solve the FULL nonlinear single-wavelength problem

  Div S = 0,   S = F - p cof(F)  (mu=1),   det F = 1,   F = F0 + grad w,
  base X2=H clamped,  top X2=0 traction-free + nonlinear Young-Laplace,

pin the fundamental modal amplitude A (lambda floats), and trace lambda(A).
Near onset  lambda(A) = lam_c + LAM2 * A^2 + ... ;  eps = 1-lambda, so
  eps - eps_c = -LAM2 A^2.   LAM2>0 => SUBCRITICAL (crease-like);  LAM2<0 => super.

Fourier collocation in X1 (one wavelength), Chebyshev in X2.  Newton via scipy.
"""
import numpy as np
import weakly_nonlinear as wn


def _fd_jac(res, U, A, r0, h=1e-7):
    J = np.empty((len(r0), len(U)))
    for j in range(len(U)):
        Up = U.copy(); Up[j] += h
        J[:, j] = (res(Up, A) - r0) / h
    return J


def newton(res, U0, A, iters=25, tol=1e-11):
    U = U0.copy()
    nr = np.linalg.norm(res(U, A))
    for it in range(iters):
        r = res(U, A); nr = np.linalg.norm(r)
        if nr < tol:
            break
        J = _fd_jac(res, U, A, r)
        dU = np.linalg.lstsq(J, -r, rcond=None)[0]
        s = 1.0
        for _ in range(10):                     # backtracking line search
            if np.linalg.norm(res(U + s * dU, A)) < nr:
                break
            s *= 0.5
        U = U + s * dU
    return U, np.linalg.norm(res(U, A))


def make_solver(N, M, H, k0, gam):
    X2, D, D2 = wn.operators(N, H)
    n = N + 1
    L1 = 2 * np.pi / k0
    X1 = L1 * np.arange(M) / M
    kx = np.fft.fftfreq(M, d=1.0 / M) * k0          # X1 wavenumbers

    def dX1(F):
        return np.real(np.fft.ifft(1j * kx[None, :] * np.fft.fft(F, axis=1), axis=1))

    def unpack(U):
        a = U[:3 * n * M]
        w1 = a[0:n * M].reshape(n, M)
        w2 = a[n * M:2 * n * M].reshape(n, M)
        p = a[2 * n * M:3 * n * M].reshape(n, M)
        lam = U[3 * n * M]
        return w1, w2, p, lam

    def residual(U, A):
        w1, w2, p, lam = unpack(U)
        H00 = dX1(w1); H01 = D @ w1
        H10 = dX1(w2); H11 = D @ w2
        F00 = lam + H00; F01 = H01; F10 = H10; F11 = 1.0 / lam + H11
        detF = F00 * F11 - F01 * F10
        # cof(F) = [[F11,-F10],[-F01,F00]]
        S00 = F00 - p * F11
        S01 = F01 - p * (-F10)
        S10 = F10 - p * (-F01)
        S11 = F11 - p * F00
        R1 = dX1(S00) + D @ S01            # (Div S)_x
        R2 = dX1(S10) + D @ S11            # (Div S)_y
        R3 = detF - 1.0                    # incompressibility
        # --- BCs ---
        # base X2=H (row N): clamp
        R1[N, :] = w1[N, :]
        R2[N, :] = w2[N, :]
        # top X2=0 (row 0): traction-free + nonlinear Young-Laplace
        tx = lam + H00[0, :]; ty = H10[0, :]      # deformed-surface tangent
        nrm = np.sqrt(tx**2 + ty**2)
        taux = tx / nrm; tauy = ty / nrm
        tYLx = -gam * dX1(taux[None, :])[0]       # -gamma d(tau)/dX1  (sign set by linear limit)
        tYLy = -gam * dX1(tauy[None, :])[0]
        R1[0, :] = S01[0, :] - tYLx               # shear traction = YL_x
        R2[0, :] = S11[0, :] - tYLy               # normal traction = YL_y
        # pin: fundamental cos(k0 X1) component of surface w2 = A
        pin = (2.0 / M) * np.sum(w2[0, :] * np.cos(k0 * X1)) - A
        return np.concatenate([R1.ravel(), R2.ravel(), R3.ravel(), [pin]])

    return residual, (X2, D, X1, n)


def branch(N, M, H, k0, gam, lam_c, f1, f2, f3, amps):
    res, (X2, D, X1, n) = make_solver(N, M, H, k0, gam)
    pi0 = 1.0 / lam_c**2
    lams = []
    Uprev = None
    for A in amps:
        # initial guess: trivial + A*eigenmode
        w1 = A * np.outer(f1, np.sin(k0 * X1))
        w2 = A * np.outer(f2, np.cos(k0 * X1))
        p = pi0 + A * np.outer(f3, np.cos(k0 * X1))
        U0 = np.concatenate([w1.ravel(), w2.ravel(), p.ravel(), [lam_c]])
        if Uprev is not None:
            U0 = Uprev.copy(); U0[-1] = lams[-1]
        U, rn = newton(res, U0, A)
        lam = U[3 * n * M]
        lams.append(lam)
        Uprev = U
        print(f"    A={A:.4f}: lam={lam:.6f}  eps={1-lam:.6f}  |res|={rn:.1e}", flush=True)
    return np.array(lams)


if __name__ == "__main__":
    N, M, H = 16, 10, 4.0
    gam = 2.0
    X2, D, D2 = wn.operators(N, H)
    # selected k0 and lam_c from validated linear solver
    best = None
    for k in np.linspace(0.3, 1.6, 20):
        lc, sv = wn.critical_lambda(N, D, D2, k, gam)
        if lc is not None and (best is None or (1 - lc) < best[0]):
            best = (1 - lc, k, lc)
    eps_c, k0, lam_c = best
    f1, f2, f3 = wn.eigenmode(N, D, D2, lam_c, k0, gam)
    f1, f2, f3 = f1 / f2[0], f2 / f2[0], f3 / f2[0]
    print(f"gam={gam} (gbar={gam/H:.2f}): k0={k0:.3f}, lam_c={lam_c:.5f}, eps_c={eps_c:.4f}")
    print("  tracing branch lambda(A):")
    amps = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
    lams = branch(N, M, H, k0, gam, lam_c, f1, f2, f3, amps)
    # fit lam = lam_c + LAM2 A^2 (use small-A points)
    c = np.polyfit(amps**2, lams, 1)
    LAM2 = c[0]; lam0 = c[1]
    print(f"  fit: lam(A) = {lam0:.5f} + ({LAM2:+.3f}) A^2   (linear lam_c={lam_c:.5f})")
    print(f"  => eps - eps_c = {-LAM2:+.3f} A^2   [{'SUBcritical (crease)' if LAM2>0 else 'SUPERcritical (wrinkle)'}]")
