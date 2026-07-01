#!/usr/bin/env python3
"""Extract the critical NOMINAL field and wavelength from the EM unstrained sweep
(em_g*.io0.exo).  Applied top voltage Phi(t)=Vmax*(t-N)/RAMP (clamped to [0,Vmax]),
nominal field Etilde=Phi/H_f, normalized Etilde*sqrt(eps/mu) with eps=mu=1.
Onset: top-surface amplitude (RMS of the cubic-detrended profile /H) first exceeds
ONSET_THR above the noise floor (sharp undamped bifurcation).  Sides are confined
(u_x=0), so the in-plane stretch stays 1 and the current wavelength = reference."""
import glob, json, os, re
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
H = 4.0
VMAX = 18.0; N = 50.0; RAMP = 1000.0
ONSET_THR = 0.03        # EM wrinkle amplitudes are small; displacement is in the io1 file


def analyze(f):
    with Dataset(f) as ds:
        cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
        DY = np.array(ds["vals_nod_var2"][:]); t = np.array(ds["time_whole"][:])
    top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
    o = np.argsort(cx[top]); top = top[o]; xr = cx[top]
    L = xr.max() - xr.min(); I = (xr > 0.08 * L) & (xr < 0.92 * L)
    Phi = np.clip((t - N) / RAMP, 0.0, 1.0) * VMAX
    Etil = Phi / H                                      # nominal field (eps=mu=1)
    A = np.zeros(len(t)); prof = []
    for k in range(len(t)):
        y = cy[top] + DY[k, top]
        yd = y - np.polyval(np.polyfit(xr, y, 3), xr)
        A[k] = np.sqrt(np.mean(yd[I] ** 2)) / H; prof.append(yd)
    on = np.where((A > ONSET_THR) & (Etil > 0.1))[0]
    if not len(on):
        return None, np.nan, A.max()
    k0 = on[0]
    # wavelength (FFT) at onset; sides confined -> current = reference
    xi = xr[I]; dx = xi[1] - xi[0]; lams = []
    for k in range(k0, min(k0 + 4, len(t))):
        yd = prof[k][I]; yd = yd - yd.mean()
        F = np.abs(np.fft.rfft(yd, n=8 * len(yd))); fr = np.fft.rfftfreq(8 * len(yd), d=dx)
        j = 1 + int(np.argmax(F[1:]))
        if fr[j] > 0: lams.append(1.0 / fr[j])
    lH = (np.median(lams) / H) if lams else np.nan
    return float(Etil[k0]), lH, A.max()


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "em_g*.io1.exo"))):   # io1 = displacement (D_X,D_Y)
        gbar = int(re.search(r"em_g(\d+)", f).group(1)) / 100.0
        Ec, lH, amax = analyze(f)
        rows.append({"gbar": gbar, "Ec": round(Ec, 3) if Ec else None,
                     "lH": round(float(lH), 3) if not np.isnan(lH) else None,
                     "amp_max": round(float(amax), 4)})
        tag = f"Ec(nominal)={Ec:.3f} l/H={lH:.2f}" if Ec else f"no onset (amax={amax:.3f})"
        print(f"  gbar={gbar:>5}: {tag}")
    json.dump(rows, open(os.path.join(HERE, "em_onset.json"), "w"), indent=2)
    print(f"wrote em_onset.json ({len(rows)} runs)")


if __name__ == "__main__":
    main()
