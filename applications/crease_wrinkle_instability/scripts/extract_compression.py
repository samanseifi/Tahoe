#!/usr/bin/env python3
"""
Extract the wrinkling onset strain eps_c and wavelength from a mechanical
lateral-compression run (compress_single_2D.xml / compress_bilayer_2D.xml).

Nominal compressive strain vs simulation time: eps(t) = rate * t, where
rate = (uL + uR)/Lx0 per unit time = (|vL| + |vR|)/Lx0  (default grips
+-0.08 on Lx0=40  ->  rate = 0.16/40 = 0.004).

Onset eps_c = strain at which the top-surface peak-to-peak D_Y first exceeds
a tolerance (above numerical noise). Wavelength from FFT of the top-edge
D_Y profile a few frames past onset.

Usage:
  python3 extract_compression.py --exo compress_single.io1.exo \
        --ytop 4.0 --rate 0.004 --Lx 40 --tol 0.02
"""
import argparse
import numpy as np
from netCDF4 import Dataset


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exo", required=True)
    ap.add_argument("--ytop", type=float, default=4.0)
    ap.add_argument("--rate", type=float, default=0.004, help="d(eps)/d(t)")
    ap.add_argument("--Lx", type=float, default=40.0)
    ap.add_argument("--tol", type=float, default=0.02, help="onset amplitude (abs)")
    args = ap.parse_args()

    with Dataset(args.exo, "r") as ds:
        cx = np.asarray(ds.variables["coordx"][:])
        cy = np.asarray(ds.variables["coordy"][:])
        t = np.asarray(ds.variables["time_whole"][:])
        DX = np.asarray(ds.variables["vals_nod_var1"][:])
        DY = np.asarray(ds.variables["vals_nod_var2"][:])

    top = np.where(np.isclose(cy, args.ytop, atol=1e-6))[0]
    order = np.argsort(cx[top]); top = top[order]
    x0 = cx[top]
    amp = DY[:, top].max(axis=1) - DY[:, top].min(axis=1)
    eps = args.rate * t

    # onset: first frame above tol
    k = int(np.argmax(amp > args.tol))
    if amp[k] <= args.tol:
        print(f"no onset reached (max amp {amp.max():.2e} < tol {args.tol}); "
              f"max eps={eps[-1]:.3f}")
        return
    eps_c = eps[k]
    print(f"onset: frame {k}, eps_c = {eps_c:.4f}  (t={t[k]:.2f}, amp={amp[k]:.3e})")

    # wavelength: a few frames past onset (clear signal, pre-inversion)
    kf = min(k + 3, len(t) - 1)
    prof = DY[kf, top]
    prof = prof - prof.mean()
    # interior window to avoid grip ends
    xi = x0[(x0 > 0.1 * args.Lx) & (x0 < 0.9 * args.Lx)]
    pi = prof[(x0 > 0.1 * args.Lx) & (x0 < 0.9 * args.Lx)]
    n = len(pi)
    dx = (xi[-1] - xi[0]) / (n - 1)
    fft = np.abs(np.fft.rfft(pi - pi.mean()))
    freqs = np.fft.rfftfreq(n, d=dx)
    j = 1 + int(np.argmax(fft[1:]))
    lam = 1.0 / freqs[j] if freqs[j] > 0 else float("nan")
    print(f"wavelength at eps={eps[kf]:.3f}: l_c = {lam:.2f}  (l_c/H = {lam/args.ytop:.2f})")
    print(f"  [theory case1 gamma_bar=0.5: eps_c~0.594, l_c/H~4.1]")


if __name__ == "__main__":
    main()
