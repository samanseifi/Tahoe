#!/usr/bin/env python3
"""Extract electromechanical wrinkle onset (critical nominal field Etilde_c and
wavelength) from the stage-2 explicit voltage-ramp runs.

Voltage schedule 2: Phi = 20*s(t), s=0 for t<=80 then linear to 1 at t=680, so
the normalized nominal field is  Etilde*sqrt(eps/mu) = Phi/H = 5*s(t),  H=4.
Onset = first frame where the detrended top-surface undulation exceeds A_ON."""
import os
import numpy as np
from netCDF4 import Dataset

ROOT = os.path.join(os.path.dirname(__file__), "..")
H = 4.0
A_ON = 0.30          # peak-to-peak onset threshold (above explicit-noise floor)


def Efield(t):
    s = np.clip((t - 80.0) / 600.0, 0.0, 1.0)
    return 5.0 * s


def analyze(path):
    with Dataset(path) as ds:
        cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
        t = np.array(ds["time_whole"][:])
        DX = np.array(ds["vals_nod_var1"][:]); DY = np.array(ds["vals_nod_var2"][:])
    top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
    x0 = cx[top]
    Lx = x0.max() - x0.min()
    interior = (x0 > x0.min() + 0.12 * Lx) & (x0 < x0.max() - 0.12 * Lx)
    amp = np.zeros(len(t))
    for k in range(len(t)):
        xc = x0[interior] + DX[k, top][interior]
        yc = cy[top][interior] + DY[k, top][interior]
        srt = np.argsort(xc); xc = xc[srt]; yc = yc[srt]
        yd = yc - np.polyval(np.polyfit(xc, yc, 2), xc)   # remove macroscopic bending
        amp[k] = np.ptp(yd)
    E = Efield(t)
    # onset: first frame with field on (E>0.05) and amp>A_ON
    on = np.where((amp > A_ON) & (E > 0.05))[0]
    if not len(on):
        return None
    k = on[0]
    Ec = E[k]
    # wavelength at a developed frame (amp in [A_ON, 5x A_ON]); zero-padded FFT
    dev = np.where((amp > A_ON) & (amp < 8 * A_ON) & (E > 0.05))[0]
    kk = dev[-1] if len(dev) else k
    xc = x0[interior] + DX[kk, top][interior]
    yc = cy[top][interior] + DY[kk, top][interior]
    srt = np.argsort(xc); xc = xc[srt]; yc = yc[srt]
    yd = yc - np.polyval(np.polyfit(xc, yc, 2), xc)
    n = len(yd); dx = (xc[-1] - xc[0]) / (n - 1)
    F = np.abs(np.fft.rfft(yd, n=4 * n)); fr = np.fft.rfftfreq(4 * n, dx)
    j = 1 + int(np.argmax(F[1:]))
    lam = (1.0 / fr[j]) if fr[j] > 0 else np.nan
    return Ec, lam / H, t[k]


CASES = []
for gb, g in [(0.5, 2), (1.0, 4), (2.0, 8), (5.0, 20)]:
    CASES.append(("compression", 0.9, gb, f"em_compression/emc_s2_comp_g{g}.io1.exo"))
    CASES.append(("stretch",     1.2, gb, f"em_stretch/emc_s2_stretch_g{g}.io1.exo"))

print(f"{'case':>12} {'lam':>4} {'gbar':>5} {'Ec(nominal)':>11} {'l/H':>6} {'t_onset':>8}")
results = {}
for kind, lam, gb, rel in CASES:
    path = os.path.join(ROOT, rel)
    if not os.path.exists(path):
        print(f"{kind:>12} {lam:>4} {gb:>5} {'(no exo)':>11}"); continue
    r = analyze(path)
    if r is None:
        print(f"{kind:>12} {lam:>4} {gb:>5} {'(no onset)':>11}"); continue
    Ec, lH, ton = r
    results[(kind, gb)] = (Ec, lH)
    print(f"{kind:>12} {lam:>4} {gb:>5} {Ec:>11.3f} {lH:>6.2f} {ton:>8.1f}")

import json
json.dump({f"{k[0]}_{k[1]}": v for k, v in results.items()},
          open(os.path.join(ROOT, "em_fe_onset.json"), "w"), indent=2)
print("\nwrote em_fe_onset.json")
