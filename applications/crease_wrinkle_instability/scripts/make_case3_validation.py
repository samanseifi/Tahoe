#!/usr/bin/env python3
"""
Case 3 FE validation figure: overlay mechanical-compression FE onset strains
and wavelengths on the bilayer linear-stability prediction, for the FE
geometry (Lx=40, H1=8 substrate, H2=1 film). Produces pics/case3_fe.pdf.

FE runs (multilayer_experiments/):
  compress_bilayer_2D.io1.exo : mu2/mu1=10, gamma_int=0
  cb_r5.io1.exo               : mu2/mu1=5,  gamma_int=0
  cb_r20.io1.exo              : mu2/mu1=20, gamma_int=0
  cb_gi5.io1.exo              : mu2/mu1=10, gamma_int=5  (elastocapillary shift)
"""
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
import nlayer_stability as nl
sys.path.insert(0, os.path.dirname(__file__))
import plotstyle as ps
ps.apply_style()

EXP = os.path.join(os.path.dirname(__file__), "..", "fe_tests", "5_mechanical_multilayer")
PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")
PAIR = ps.PAIR; REF = "#8a8f99"


def fe_onset(exo, ytop=9.0, Lx=40.0, rate=0.002, A_on=0.1):
    """Onset = strain at which the detrended top-surface undulation first
    reaches A_on (10%% of the film thickness H2=1), a clean monotonic
    criterion; developed wavelength taken just past onset. Returns (None,nan)
    if onset is not reached within the run."""
    with Dataset(exo) as ds:
        cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
        t = np.array(ds["time_whole"][:]); DY = np.array(ds["vals_nod_var2"][:])
    row = np.where(np.isclose(cy, ytop, atol=1e-6))[0]; row = row[np.argsort(cx[row])]
    x = cx[row]; m = (x > 0.1 * Lx) & (x < 0.9 * Lx); xi = x[m]
    eps = rate * t
    amp = np.array([np.ptp(p - np.polyval(np.polyfit(xi, p, 2), xi))
                    for p in (DY[k, row][m] for k in range(len(t)))])
    above = np.where((amp > A_on) & (eps > 0.04))[0]
    if len(above) == 0:
        return None, np.nan
    k = int(above[0]); ec = eps[k]
    kd = int(np.argmin(np.abs(eps - (ec + 0.02))))
    p = DY[kd, row][m]; p = p - np.polyval(np.polyfit(xi, p, 2), xi)
    n = len(p); dx = (xi[-1] - xi[0]) / (n - 1)
    F = np.abs(np.fft.rfft(p)); fr = np.fft.rfftfreq(n, dx); j = 1 + int(np.argmax(F[1:]))
    lam = 1.0 / fr[j] if fr[j] > 0 else np.nan
    return ec, lam


def main():
    H1, H2 = 8.0, 1.0
    runs = {  # (exo, mu2, gamma_int)
        5:  ("cb_r5.io1.exo", 5.0, 0.0),
        10: ("compress_bilayer_2D.io1.exo", 10.0, 0.0),
        20: ("cb_r20.io1.exo", 20.0, 0.0),
    }
    fe_r, fe_ec, fe_lam = [], [], []
    for r in sorted(runs):
        exo, mu2, gi = runs[r]
        path = os.path.join(EXP, exo)
        if not os.path.exists(path):
            print(f"missing {exo}"); continue
        ec, lam = fe_onset(path)
        if ec is None:
            print(f"FE ratio={r}: onset not reached in run (skipped)"); continue
        fe_r.append(r); fe_ec.append(ec); fe_lam.append(lam)
        print(f"FE ratio={r}: eps_c={ec:.4f}, lambda={lam:.2f}")

    # theory curves (gamma=0)
    ratios = np.geomspace(3, 30, 40); Ks = np.linspace(0.05, 3.0, 220)
    th_ec, th_lam = [], []
    for r in ratios:
        res = nl.critical([1.0, r], [H1, H2], [0.0], 0.0, Ks)
        th_ec.append(res[0]); th_lam.append(2 * np.pi / res[1])

    # interfacial-energy point (ratio 10): FE gi=0 vs gi=5 and theory
    gi_path = os.path.join(EXP, "cb_gi5.io1.exo")
    gi_fe = fe_onset(gi_path)[0] if os.path.exists(gi_path) else None
    th_gi0 = nl.critical([1.0, 10.0], [H1, H2], [0.0], 0.0, Ks)[0]
    th_gi5 = nl.critical([1.0, 10.0], [H1, H2], [5.0], 0.0, Ks)[0]

    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.8))
    ax[0].plot(ratios, th_ec, "-", color=PAIR[0], label="linear theory")
    ps.fe_marker(ax[0], fe_r, fe_ec, color=PAIR[1], marker="o", label="FE (compression)")
    ax[0].set_ylabel(r"onset strain $\varepsilon_c$")
    ax[1].plot(ratios, th_lam, "-", color=PAIR[0], label="linear theory")
    ps.fe_marker(ax[1], fe_r, fe_lam, color=PAIR[1], marker="o", label="FE (compression)")
    ax[1].set_ylabel(r"wavelength $\ell_c$")
    for a in ax:
        a.set_xlabel(r"modulus ratio $\mu_2/\mu_1$"); a.set_xscale("log")
        a.legend(frameon=False)
    ax[0].set_title(r"Onset strain (film/substrate, $H_1{=}8,H_2{=}1$)")
    ax[1].set_title(r"Selected wavelength")
    fig.tight_layout(); fig.savefig(f"{PICS}/case3_fe.pdf"); plt.close(fig)

    print(f"\nInterfacial-energy shift (ratio 10): "
          f"theory eps_c {th_gi0:.3f} -> {th_gi5:.3f}; "
          f"FE {fe_ec[fe_r.index(10)]:.3f} -> {gi_fe if gi_fe else 'NA'}")
    print(f"wrote {PICS}/case3_fe.pdf")


if __name__ == "__main__":
    main()
