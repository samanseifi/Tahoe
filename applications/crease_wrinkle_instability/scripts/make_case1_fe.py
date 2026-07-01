#!/usr/bin/env python3
"""
Overlay the single-layer mechanical-compression FE results on the Case-1 sliding-
base theory: critical strain (ge_mech.pdf) and wrinkle wavelength (gl_mech.pdf),
with classic open markers (plotstyle.fe_marker).

Onset strains are extrapolated to zero amplitude from the symmetric-compression
sweep (roller base, scales L=80/H=4 and L=160/H=8 agree).  Wavelengths are from
the L=80 sweep (two estimators); L/H=20 best suppresses the competing domain mode,
so these are the most reliable -- the longer L=180 strip gives that sliding mode
more room and degrades the wavelength, while the onset is unchanged.  The wavelength
agrees with theory only to ~25% (see text); the onset to ~6%.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import nlayer_stability as nl
import plotstyle as ps
ps.apply_style()
import matplotlib.pyplot as plt

PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")
SOLO = ps.PALETTE[0]
REF = "#8a8f99"
FE_C = ps.PAIR[1]
H = 4.0

# FE results -- undamped symmetric-compression sweep, roller base, L=180/H=4
# (no damping -> the wrinkle develops in the bulk, not pinned at the boundary).
FE_GBAR = np.array([0.5, 1.0, 2.0, 3.5, 5.0])
FE_EPSC = np.array([0.619, 0.667, 0.720, 0.760, 0.788])   # onset (threshold above noise)
# wavelength: CURRENT (deformed) wavelength l_c/H_f = (1-eps_c) * l_ref/H_f -- the
# physically measured crest-to-crest spacing in the compressed film.  gbar<=3.5,
# where the wave fits >=3 periods; at gbar=5 the finite domain quantizes to 2 waves.
FE_LGBAR = np.array([0.5, 1.0, 2.0, 3.5])
FE_EPSC_L = np.array([0.619, 0.667, 0.720, 0.760])        # onset strain at these gbar
FE_LH = np.array([4.98, 6.37, 10.15, 14.42]) * (1.0 - FE_EPSC_L)   # -> current l/H


def main():
    # high-resolution, smooth theory curves.  LOG-spaced K resolves the small-K
    # (long-wavelength) region at high gbar, so l_c/H=2pi/K_c comes out smooth
    # (linear K spacing quantizes K_c and makes the wavelength curve jagged).
    Ks = np.geomspace(0.03, 6.0, 70)
    gb_fine = np.linspace(0.0, 5.0, 34)
    th_ec = np.array([nl.critical([1.0], [H], [], gb * H, Ks, base="sliding")[0] for gb in gb_fine])

    # (a) critical strain
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.plot(gb_fine, th_ec, "-", color=SOLO, label="linear theory", zorder=2)
    ps.fe_marker(ax, FE_GBAR, FE_EPSC, color=FE_C, marker="o", label="FE (compression)")
    ax.axhline(0.456, ls=":", c=REF, zorder=1)
    ax.text(3.6, 0.47, r"Biot ($\bar\gamma=0$)", color=REF, fontsize=10)
    ax.set_xlabel(r"elastocapillary number $\bar\gamma=\gamma/(\mu H_f)$")
    ax.set_ylabel(r"critical strain $\varepsilon_c$")
    ax.set_xlim(-0.15, 5.2)                     # margin so gbar=5 marker is not on the spine
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{PICS}/ge_mech.pdf"); plt.close(fig)

    # (b) wrinkle wavelength  (FFT estimates only) -- CURRENT (deformed) wavelength
    #     l_c/H = (2*pi/K_c) * lambda_c / H = (1-eps_c) * l_ref/H
    gb_l = np.linspace(0.1, 5.0, 34)
    th_l = []
    for gb in gb_l:
        r = nl.critical([1.0], [H], [], gb * H, Ks, base="sliding")
        th_l.append(2 * np.pi / r[1] / H * r[2])          # r[2] = lambda_c = 1 - eps_c
    th_l = np.array(th_l)
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.plot(gb_l, th_l, "-", color=SOLO, label="linear theory", zorder=2)
    ps.fe_marker(ax, FE_LGBAR, FE_LH, color=FE_C, marker="o", label="FE (compression)")
    ax.set_xlabel(r"elastocapillary number $\bar\gamma=\gamma/(\mu H_f)$")
    ax.set_ylabel(r"critical wavelength $\ell_c/H_f$")
    ax.set_xlim(-0.15, 5.2)
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{PICS}/gl_mech.pdf"); plt.close(fig)
    print("ge_mech + gl_mech (theory + %d FE markers each) done" % len(FE_GBAR))


if __name__ == "__main__":
    main()
