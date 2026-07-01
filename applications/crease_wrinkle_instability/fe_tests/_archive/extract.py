#!/usr/bin/env python3
"""Onset extraction for the clean pre-strain sweep (stage2_l*.io1.exo).
Nominal field ramp:  Etilde = Phi/H = 10 * clip((t-40)/1200, 0, 1)."""
import glob
import numpy as np
from netCDF4 import Dataset


def Efield(t):
    return np.clip((t - 40.0) / 1200.0, 0.0, 1.0) * 10.0


def amp_series(f):
    with Dataset(f) as ds:
        cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
        t = np.array(ds["time_whole"][:])
        DX = np.array(ds["vals_nod_var1"][:]); DY = np.array(ds["vals_nod_var2"][:])
    top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
    x0 = cx[top]; Lx = x0.max() - x0.min()
    I = (x0 > x0.min() + 0.12 * Lx) & (x0 < x0.max() - 0.12 * Lx)
    amp = np.zeros(len(t))
    for k in range(len(t)):
        xc = x0[I] + DX[k, top][I]; yc = cy[top][I] + DY[k, top][I]
        s = np.argsort(xc); xc = xc[s]; yc = yc[s]
        yd = yc - np.polyval(np.polyfit(xc, yc, 2), xc)
        amp[k] = np.ptp(yd)
    return t, amp, Efield(t), (x0, I, top, cx, cy, DX, DY)


def wavelength(f, k):
    with Dataset(f) as ds:
        cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
        DX = np.array(ds["vals_nod_var1"][:]); DY = np.array(ds["vals_nod_var2"][:])
    top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
    x0 = cx[top]; Lx = x0.max() - x0.min()
    I = (x0 > x0.min() + 0.12 * Lx) & (x0 < x0.max() - 0.12 * Lx)
    xc = x0[I] + DX[k, top][I]; yc = cy[top][I] + DY[k, top][I]
    s = np.argsort(xc); xc = xc[s]; yc = yc[s]
    yd = yc - np.polyval(np.polyfit(xc, yc, 2), xc)
    n = len(yd); dx = (xc[-1] - xc[0]) / (n - 1)
    F = np.abs(np.fft.rfft(yd, n=4 * n)); fr = np.fft.rfftfreq(4 * n, dx)
    j = 1 + int(np.argmax(F[1:]))
    return (1.0 / fr[j]) / 4.0 if fr[j] > 0 else np.nan


print(f"{'lam':>6} {'eps_pre':>8} {'frames':>7} {'E_max':>6} {'max_amp':>8} {'Etilde_c':>9} {'l/H':>5}")
rows = []
for f in sorted(glob.glob("stage2_l*.io1.exo")):
    lam = float(f.split("_l")[1][:3]) / 100.0
    t, amp, E, _ = amp_series(f)
    on = np.where((amp > 0.3) & (E > 0.05))[0]
    Eon = E[on[0]] if len(on) else np.nan
    lH = wavelength(f, on[0]) if len(on) else np.nan
    rows.append((lam, Eon, lH))
    print(f"{lam:>6.2f} {1-lam:>+8.2f} {len(t):>7} {E.max():>6.2f} {amp.max():>8.2f} {Eon:>9.3f} {lH:>5.2f}")
import json
json.dump([{"lam": l, "Ec": (None if np.isnan(e) else e), "lH": (None if np.isnan(w) else w)}
           for l, e, w in rows], open("sweep_onset.json", "w"), indent=2)
