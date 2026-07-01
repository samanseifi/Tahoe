#!/usr/bin/env python3
"""
Case 2 (electromechanical) linear-stability figures -> paper_draft/pics/.
Wrinkling criterion (field-free biharmonic bulk; the field couples through the
top Maxwell traction) for pre-strained films, in elec_stability.py.
JAM-style figures (plotstyle.py).

  e2_elec.pdf         : neutral-stability curves Etilde*sqrt(eps/mu) vs l/H_f at
                        several elastocapillary numbers, pre-COMPRESSION lam=0.9.
  e2_elec_stretch.pdf : the same family at pre-STRETCH lam=1.2.
  wrinkle_elec1.pdf   : critical nominal field vs gamma_bar, compression+stretch.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import plotstyle as ps
ps.apply_style()
import matplotlib.pyplot as plt
import elec_stability as es

PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")
H = 4.0


def neutral_panel(lam, gbars, fname, title):
    Ks = np.linspace(0.10, 2.4, 70)                  # higher resolution
    fig, ax = plt.subplots(figsize=(5.0, 3.9))
    for i, gbar in enumerate(gbars):
        ll, EE = es.neutral_curve(lam, gbar, H, Ks)
        if len(ll):
            o = np.argsort(ll)
            ax.plot(ll[o], EE[o], color=ps.PALETTE[i], label=fr"$\bar\gamma={gbar:g}$")
    ax.set_xlabel(r"normalized wavelength $\ell/H_f$")
    ax.set_ylabel(r"nominal field $\tilde E\,\sqrt{\epsilon/\mu}$")
    ax.set_xlim(0, 12)
    ax.legend(ncol=2, frameon=False)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(f"{PICS}/{fname}")
    plt.close(fig)
    print(f"{fname} done")


def main():
    gbars = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0]
    neutral_panel(0.8, gbars, "e2_elec.pdf",
                  r"Electromechanical wrinkling, $\lambda_{pre}=0.8$")
    neutral_panel(1.3, gbars, "e2_elec_stretch.pdf",
                  r"Electromechanical wrinkling, $\lambda_{pre}=1.3$")

    # critical nominal field vs gamma_bar (to 20): compression vs stretch
    Ks = np.linspace(0.05, 2.6, 80)
    gbg = np.linspace(0.0, 20.0, 41)
    out = {}
    for lam in (0.80, 1.30):
        out[lam] = np.array([(es.critical_field(lam, gb, H, Ks) or [np.nan])[0]
                             for gb in gbg])
        print(f"lam={lam}: Ec(0)={out[lam][0]:.3f}, Ec(20)={out[lam][-1]:.3f}")
    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    ax.plot(gbg, out[0.80], color=ps.PAIR[0],
            label=r"pre-compression ($\lambda_{pre}=0.8$)")
    ax.plot(gbg, out[1.30], color=ps.PAIR[1], ls=(0, (5, 4)),
            label=r"pre-stretch ($\lambda_{pre}=1.3$)")
    ax.set_xlabel(r"elastocapillary number $\bar\gamma=\gamma/(\mu H_f)$")
    ax.set_ylabel(r"critical nominal field $\tilde E_c\,\sqrt{\epsilon/\mu}$")
    ax.set_xlim(0, 20)
    ax.legend(frameon=False)
    ax.set_title(r"Critical field vs surface tension")
    fig.tight_layout()
    fig.savefig(f"{PICS}/wrinkle_elec1.pdf")
    plt.close(fig)
    print("wrinkle_elec1 done")


if __name__ == "__main__":
    main()
