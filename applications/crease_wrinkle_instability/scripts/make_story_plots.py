#!/usr/bin/env python3
"""
Linear-stability figures for the surface-instability story, written as PDFs
into paper_draft/pics/. Smooth curves (refined K_c, high resolution), LINES
ONLY -- markers are intentionally reserved for overlaying FEA results later.

Cases (mechanical compression; validated):
  CASE 1  single layer + surface tension          -> ew_mech, ge_mech, gl_mech
  CASE 3  thin film / thick substrate + gamma_i    -> bilayer_map
  CASE 5  N-layer mode selection (gamma_t)         -> nlayer_mode

Electromechanical cases (2,4,6) are intentionally NOT plotted (need the
Maxwell-stress interface jump + Wang-Zhao validation).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.dirname(__file__))
import nlayer_stability as nl
import plotstyle as ps
ps.apply_style()

PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")
os.makedirs(PICS, exist_ok=True)
PALETTE = ps.PALETTE
PAIR = ps.PAIR                 # paired comparisons
SOLO = ps.PALETTE[0]
REF = "#8a8f99"                # neutral grey for reference lines


def eps_at_K(mus, Hs, g_int, g_top, K, base="fixed"):
    lam = nl.crit_lambda_at_K(mus, Hs, g_int, g_top, K, base=base)
    return None if lam is None else 1.0 - lam


# ----------------------------------------------------------------------
# CASE 1 — single layer + surface tension
# ----------------------------------------------------------------------
def case1():
    mu, H = 1.0, 4.0
    # log-spaced K so long wavelengths (small K) are well sampled out to l/H=20
    Ks = np.geomspace(0.04, 6.0, 300)

    # (a) neutral-stability curves: strain vs normalized wavelength
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for gb in [0.0, 0.5, 1.0, 2.0, 3.5, 5.0]:
        g = gb * mu * H
        ll, ee = [], []
        for K in Ks:
            e = eps_at_K([mu], [H], [], g, K, base="sliding")
            if e is not None:
                ll.append(2 * np.pi / K / H * (1 - e)); ee.append(e)   # CURRENT (deformed) wavelength
        ax.plot(ll, ee, label=fr"$\bar\gamma={gb}$")
    ax.axhline(0.456, ls=":", c=REF, lw=1.4)
    ax.text(0.4, 0.467, r"Biot $0.456$", color=REF, fontsize=10, va="bottom", ha="left")
    ax.set_xlabel(r"normalized wavelength $\ell/H_f$")
    ax.set_ylabel(r"compressive strain $\varepsilon$")
    ax.set_xlim(0, 20); ax.set_ylim(0.43, 0.92)
    ax.legend(ncol=2, frameon=False, loc="upper right")
    ax.set_title(r"Single layer $+$ surface tension")
    fig.tight_layout(); fig.savefig(f"{PICS}/ew_mech.pdf"); plt.close(fig)

    # (b,c) critical strain and wavelength vs elastocapillary number
    gb_fine = np.linspace(0.0, 5.0, 51)
    epsc, lc = [], []
    for gb in gb_fine:
        res = nl.critical([mu], [H], [], gb * mu * H, Ks, base="sliding")
        epsc.append(res[0]); lc.append(2 * np.pi / res[1] / H * res[2])   # current: * lambda_c
    for fname, y, ylab in [("ge_mech.pdf", epsc, r"critical strain $\varepsilon_c$"),
                           ("gl_mech.pdf", lc, r"critical wavelength $\ell_c/H_f$")]:
        fig, ax = plt.subplots(figsize=(4.4, 3.4))
        ax.plot(gb_fine, y, "-", color=SOLO)
        if fname == "ge_mech.pdf":
            ax.axhline(0.456, ls=":", c=REF); ax.text(3.7, 0.47, r"Biot ($\bar\gamma=0$)", color=REF, fontsize=10)
        ax.set_xlabel(r"elastocapillary number $\bar\gamma=\gamma/(\mu H_f)$")
        ax.set_ylabel(ylab)
        fig.tight_layout(); fig.savefig(f"{PICS}/{fname}"); plt.close(fig)
    print("CASE 1: eps_c(0)=%.3f  eps_c(2)=%.3f" % (epsc[0], epsc[-1]))


# ----------------------------------------------------------------------
# CASE 3 — thin film on thick substrate, with interfacial tension
# ----------------------------------------------------------------------
def case3():
    H1, H2, g_top = 20.0, 1.0, 0.5          # thick soft substrate, thin film
    ratios = np.geomspace(2.0, 50.0, 34)     # film/substrate regime
    Ks = np.linspace(0.03, 2.0, 150)
    out = {}
    for gi in (0.0, 2.0):
        ec, lc = [], []
        for r in ratios:
            res = nl.critical([1.0, r], [H1, H2], [gi], g_top, Ks)
            ec.append(res[0]); lc.append(2 * np.pi / res[1] / H2)
        out[gi] = (np.array(ec), np.array(lc))
    l_classic = 2 * np.pi * (ratios / 3.0) ** (1 / 3.0)   # /h_f, h_f=H2=1

    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.7))
    for gi, col in ((0.0, PAIR[0]), (2.0, PAIR[1])):
        lab = r"$\gamma_i=0$" if gi == 0 else r"$\gamma_i=2$"
        ax[0].plot(ratios, out[gi][0], "-", color=col, label=lab)
        ax[1].plot(ratios, out[gi][1], "-", color=col, label=lab)
    ax[1].plot(ratios, l_classic, ":", color=REF, lw=1.6,
               label=r"$2\pi(E_f/3E_s)^{1/3}$")
    ax[0].set_ylabel(r"critical strain $\varepsilon_c$")
    ax[1].set_ylabel(r"critical wavelength $\ell_c/H_2$")
    for a in ax:
        a.set_xlabel(r"modulus ratio $\mu_2/\mu_1$")
        a.set_xscale("log"); a.legend(fontsize=9, frameon=False)
    ax[0].set_title("Thin film / thick substrate: onset")
    ax[1].set_title("wavelength (vs classical law)")
    fig.tight_layout(); fig.savefig(f"{PICS}/bilayer_map.pdf"); plt.close(fig)
    print("CASE 3: ratio=10  eps_c(gi=0)=%.3f gi=2=%.3f  l/H2=%.2f"
          % (out[0.0][0][np.argmin(abs(ratios-10))],
             out[2.0][0][np.argmin(abs(ratios-10))],
             out[0.0][1][np.argmin(abs(ratios-10))]))


# ----------------------------------------------------------------------
# CASE 5 — N-layer mode selection (buried buckling vs surface)
# ----------------------------------------------------------------------
def case5():
    mus, Hs, g_int = [1.0, 20.0, 1.0], [2.0, 1.0, 2.0], [0.0, 0.0]
    Ks = np.linspace(0.05, 4.0, 150)

    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.9))
    for g_top, c in ((0.0, PAIR[0]), (30.0, PAIR[1])):
        res = nl.critical(mus, Hs, g_int, g_top, Ks)
        Y, f2, _ = nl.mode_shape(res[2], res[1], mus, Hs, g_int, g_top, npts=80)
        ax[0].plot(f2, Y, color=c, label=fr"$\gamma_t={g_top:.0f}$ ($\varepsilon_c={res[0]:.3f}$)")
    ax[0].axhspan(2.0, 3.0, color="0.85", alpha=0.6)
    ax[0].text(0.04, 2.5, "stiff sheet", fontsize=9, rotation=90, va="center")
    ax[0].set_xlabel(r"normalized eigenmode $f_2$")
    ax[0].set_ylabel(r"height $X_2$"); ax[0].legend(fontsize=9, loc="lower left", frameon=False)
    ax[0].set_title("Buckling mode shape")

    gtops = np.linspace(0, 40, 41)
    top_d, int_d = [], []
    for gt in gtops:
        res = nl.critical(mus, Hs, g_int, gt, Ks)
        _, _, bd = nl.mode_shape(res[2], res[1], mus, Hs, g_int, gt)
        int_d.append(max(bd[:-1])); top_d.append(bd[-1])
    ax[1].plot(gtops, top_d, "-", color=PAIR[0], label="top surface")
    ax[1].plot(gtops, int_d, "-", color=PAIR[1], label="buried interface (max)")
    ax[1].set_xlabel(r"top surface tension $\gamma_t$")
    ax[1].set_ylabel("normalized deflection")
    ax[1].set_title("Surface tension selects the failing interface")
    ax[1].legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(f"{PICS}/nlayer_mode.pdf"); plt.close(fig)
    print("CASE 5: top deflection %.3f -> %.3f (gamma_t 0->40)" % (top_d[0], top_d[-1]))


def demo_per_interface():
    """Show interfaces can carry DIFFERENT tensions (not forced equal)."""
    mus, Hs = [1.0, 20.0, 1.0], [2.0, 1.0, 2.0]
    Ks = np.linspace(0.05, 4.0, 200)
    print("\nPer-interface tension demo (3 layers, 2 interfaces, gamma_t=0):")
    for gi in ([0, 0], [5, 0], [0, 5], [5, 5]):
        res = nl.critical(mus, Hs, gi, 0.0, Ks)
        print(f"  gamma_int={gi}: eps_c={res[0]:.4f}  l_c={2*np.pi/res[1]:.2f}")


if __name__ == "__main__":
    case1(); case3(); case5(); demo_per_interface()
    print("\nFigures ->", os.path.normpath(PICS))
