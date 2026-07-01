#!/usr/bin/env python3
"""
Analyze a bilayer DE run: extract peak-to-peak D_Y amplitude on the TOP
surface (y = H1+H2) and on the internal INTERFACE (y = H1) versus applied
voltage, for the interface-OFF and interface-ON cases, and report the
bifurcation threshold (voltage at which the amplitude first exceeds a
small tolerance).

Usage:
    python3 analyze_bilayer.py --off bilayer_intOFF.io1.exo \
                               --on  bilayer_intON.io1.exo \
                               --H1 2 --H2 2 --vdot 0.05 --out bilayer_amp.png
"""
import argparse
import numpy as np
from netCDF4 import Dataset

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def read_row(ds, coordy, yval):
    idx = np.where(np.isclose(coordy, yval, atol=1e-6))[0]
    DX = np.asarray(ds.variables["vals_nod_var1"][:])[:, idx]
    DY = np.asarray(ds.variables["vals_nod_var2"][:])[:, idx]
    return DX, DY


def amp_history(exo, y_top, y_int):
    with Dataset(exo, "r") as ds:
        coordy = np.asarray(ds.variables["coordy"][:])
        t = np.asarray(ds.variables["time_whole"][:])
        _, DYt = read_row(ds, coordy, y_top)
        _, DYi = read_row(ds, coordy, y_int)
    amp_top = DYt.max(axis=1) - DYt.min(axis=1)
    amp_int = DYi.max(axis=1) - DYi.min(axis=1)
    return t, amp_top, amp_int


def threshold(t, amp, vdot, tol):
    k = np.argmax(amp > tol)
    if amp[k] <= tol:
        return None
    return vdot * t[k]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--off", required=True)
    ap.add_argument("--on", required=True)
    ap.add_argument("--H1", type=float, default=2.0)
    ap.add_argument("--H2", type=float, default=2.0)
    ap.add_argument("--vdot", type=float, default=0.05)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--out", default="bilayer_amp.png")
    args = ap.parse_args()

    y_top = args.H1 + args.H2
    y_int = args.H1

    to, at_top, at_int = amp_history(args.off, y_top, y_int)
    tn, an_top, an_int = amp_history(args.on, y_top, y_int)

    print(f"{'':12} {'V_crit(top)':>12} {'V_crit(int)':>12}")
    for lbl, t, atop, aint in (("interface OFF", to, at_top, at_int),
                               ("interface ON ", tn, an_top, an_int)):
        vt = threshold(t, atop, args.vdot, args.tol)
        vi = threshold(t, aint, args.vdot, args.tol)
        print(f"{lbl:12} {('%.2f'%vt) if vt else '   --':>12} "
              f"{('%.2f'%vi) if vi else '   --':>12}")

    # report amplitudes at a few matched voltages
    print("\n  V    ampTop(OFF)  ampTop(ON)   ampInt(OFF)  ampInt(ON)")
    for V in (8, 10, 12, 14, 16, 18):
        def at(t, a):
            k = int(np.argmin(np.abs(args.vdot * t - V)))
            return a[k]
        print(f"{V:4d}  {at(to,at_top):10.3e}  {at(tn,an_top):10.3e}  "
              f"{at(to,at_int):10.3e}  {at(tn,an_int):10.3e}")

    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ax[0].semilogy(args.vdot * to, at_top, "-", label="interface OFF")
        ax[0].semilogy(args.vdot * tn, an_top, "--", label="interface ON")
        ax[0].set_title("top surface (y=H1+H2)")
        ax[1].semilogy(args.vdot * to, at_int, "-", label="interface OFF")
        ax[1].semilogy(args.vdot * tn, an_int, "--", label="interface ON")
        ax[1].set_title("interface (y=H1)")
        for a in ax:
            a.set_xlabel("voltage V"); a.set_ylabel("peak-to-peak $D_Y$")
            a.grid(True, which="both", alpha=0.3); a.legend()
        fig.tight_layout()
        fig.savefig(args.out, dpi=130)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
