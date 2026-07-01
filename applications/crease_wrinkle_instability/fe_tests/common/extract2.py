#!/usr/bin/env python3
"""Onset extraction for the refined extended sweep (run_sweep2.py).
Writes two JSON files consumed by the figure scripts:
  gamma_onset.json : critical nominal field vs gbar at lambda=1  (stage2_l100_g*)
  lam_onset.json   : critical nominal field vs lambda at gbar=2  (stage2_l*_g020)
Nominal field ramp (unchanged): Etilde = Phi/H = 10*clip((t-40)/1200, 0, 1)."""
import glob
import json
import os
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
H = 4.0
THR = 0.3          # surface-amplitude onset threshold (units of length)


def Efield(t):
    return np.clip((t - 40.0) / 1200.0, 0.0, 1.0) * 10.0


def _top(ds):
    cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
    top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
    x0 = cx[top]; Lx = x0.max() - x0.min()
    I = (x0 > x0.min() + 0.12 * Lx) & (x0 < x0.max() - 0.12 * Lx)
    return cx, cy, top, x0, I


def onset(f, thr=THR):
    with Dataset(f) as ds:
        t = np.array(ds["time_whole"][:])
        DX = np.array(ds["vals_nod_var1"][:]); DY = np.array(ds["vals_nod_var2"][:])
        cx, cy, top, x0, I = _top(ds)
    amp = np.zeros(len(t))
    for k in range(len(t)):
        xc = x0[I] + DX[k, top][I]; yc = cy[top][I] + DY[k, top][I]
        s = np.argsort(xc); xc, yc = xc[s], yc[s]
        yd = yc - np.polyval(np.polyfit(xc, yc, 2), xc)
        amp[k] = np.ptp(yd)
    E = Efield(t)
    on = np.where((amp > thr) & (E > 0.05))[0]
    if not len(on):
        return np.nan, np.nan, amp.max(), len(t)
    k = on[0]
    # wavelength via FFT of the top-surface deviation at onset
    xc = x0[I] + DX[k, top][I]; yc = cy[top][I] + DY[k, top][I]
    s = np.argsort(xc); xc, yc = xc[s], yc[s]
    yd = yc - np.polyval(np.polyfit(xc, yc, 2), xc)
    n = len(yd); dx = (xc[-1] - xc[0]) / (n - 1)
    F = np.abs(np.fft.rfft(yd, n=4 * n)); fr = np.fft.rfftfreq(4 * n, dx)
    j = 1 + int(np.argmax(F[1:]))
    lH = (1.0 / fr[j]) / 4.0 / H if fr[j] > 0 else np.nan
    return E[k], lH, amp.max(), len(t)


def main():
    os.chdir(HERE)
    # gamma sweep at lambda=1
    grows = []
    print("== gamma sweep (lambda=1) ==")
    print(f"{'gbar':>5} {'frames':>6} {'amp_max':>7} {'Etilde_c':>9} {'l/H':>5}")
    for f in sorted(glob.glob("stage2_l100_g[0-9][0-9][0-9].io1.exo")):
        gtag = f.split("_g")[1][:3]
        gb = int(gtag) / 10.0
        Ec, lH, amax, nf = onset(f)
        grows.append({"gbar": gb, "Ec": None if np.isnan(Ec) else Ec,
                      "lH": None if np.isnan(lH) else lH})
        print(f"{gb:>5} {nf:>6} {amax:>7.2f} {Ec:>9.3f} {lH:>5.2f}")
    grows.sort(key=lambda r: r["gbar"])
    json.dump(grows, open("gamma_onset.json", "w"), indent=2)

    # lambda sweep at gbar=2 (include l100_g020 from the gamma sweep)
    lrows = []
    print("\n== lambda sweep (gbar=2) ==")
    print(f"{'lam':>5} {'frames':>6} {'amp_max':>7} {'Etilde_c':>9} {'E2c':>6} {'l/H':>5}")
    for f in sorted(glob.glob("stage2_l[0-9][0-9][0-9]_g020.io1.exo")):
        lam = int(f.split("_l")[1][:3]) / 100.0
        Ec, lH, amax, nf = onset(f)
        E2 = lam * Ec
        lrows.append({"lam": lam, "Ec": None if np.isnan(Ec) else Ec,
                      "E2c": None if np.isnan(E2) else E2,
                      "lH": None if np.isnan(lH) else lH})
        print(f"{lam:>5} {nf:>6} {amax:>7.2f} {Ec:>9.3f} {E2:>6.3f} {lH:>5.2f}")
    lrows.sort(key=lambda r: r["lam"])
    json.dump(lrows, open("lam_onset.json", "w"), indent=2)
    print("\nwrote gamma_onset.json, lam_onset.json")


if __name__ == "__main__":
    main()
