#!/usr/bin/env python3
"""Path-aware onset extraction for the organized fe_tests/ layout.
  gamma sweep (lambda=1): ../2_em_no_prestrain/stage2_l100_g[0-9][0-9][0-9].io1.exo
                          -> ../2_em_no_prestrain/gamma_onset.json
  pre-strain sweep (gbar=2): stage2_l*_g020.io1.exo across cases 2,3,4
                          -> ./lam_onset.json
Nominal field ramp:  Etilde = Phi/H = 10*clip((t-40)/1200, 0, 1),  H=4."""
import glob, json, os
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))      # fe_tests/common
FT = os.path.dirname(HERE)                              # fe_tests
H, THR = 4.0, 0.3


def Efield(t):
    return np.clip((t - 40.0) / 1200.0, 0.0, 1.0) * 10.0


def onset(f, thr=THR):
    with Dataset(f) as ds:
        t = np.array(ds["time_whole"][:])
        cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
        DX = np.array(ds["vals_nod_var1"][:]); DY = np.array(ds["vals_nod_var2"][:])
    top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
    x0 = cx[top]; Lx = x0.max() - x0.min()
    I = (x0 > x0.min() + 0.12 * Lx) & (x0 < x0.max() - 0.12 * Lx)
    amp = np.array([np.ptp((cy[top][I] + DY[k, top][I])
                    - np.polyval(np.polyfit(x0[I] + DX[k, top][I],
                                            cy[top][I] + DY[k, top][I], 2),
                                 x0[I] + DX[k, top][I])) for k in range(len(t))])
    E = Efield(t); on = np.where((amp > thr) & (E > 0.05))[0]
    return (E[on[0]] if len(on) else np.nan)


def main():
    # gamma sweep, lambda=1
    g = []
    for f in sorted(glob.glob(os.path.join(FT, "2_em_no_prestrain",
                                           "stage2_l100_g[0-9][0-9][0-9].io1.exo"))):
        gb = int(f.split("_g")[1][:3]) / 10.0
        Ec = onset(f); g.append({"gbar": gb, "Ec": None if np.isnan(Ec) else Ec})
        print(f"  gbar={gb:>5}: Etilde_c={Ec:.3f}")
    g.sort(key=lambda r: r["gbar"])
    json.dump(g, open(os.path.join(FT, "2_em_no_prestrain", "gamma_onset.json"), "w"), indent=2)
    # pre-strain sweep, gbar=2, across cases 2/3/4
    L = []
    for case in ("2_em_no_prestrain", "3_em_precompressed", "4_em_prestretched"):
        for f in sorted(glob.glob(os.path.join(FT, case, "stage2_l[0-9][0-9][0-9]_g020.io1.exo"))):
            lam = int(f.split("_l")[1][:3]) / 100.0
            Ec = onset(f)
            L.append({"lam": lam, "Ec": None if np.isnan(Ec) else Ec,
                      "E2c": None if np.isnan(Ec) else lam * Ec})
            print(f"  lam={lam:>4}: Etilde_c={Ec:.3f}")
    L.sort(key=lambda r: r["lam"])
    json.dump(L, open(os.path.join(HERE, "lam_onset.json"), "w"), indent=2)
    print("wrote gamma_onset.json (case 2) and lam_onset.json (common)")


if __name__ == "__main__":
    main()
