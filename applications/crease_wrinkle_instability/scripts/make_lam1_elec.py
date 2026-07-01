#!/usr/bin/env python3
"""
Replicate the free-standing electromechanical buckling at zero pre-strain
(lambda = 1), after Yang, Zhao & Sharma (JAM 2017, 84:031008).

For an unconstrained (traction-free, electroded) incompressible neo-Hookean
block held laterally at stretch lambda=1, the antisymmetric buckling critical
nominal field follows in closed form (their Eq. 72):

    [1 + (1/2)(Etilde sqrt(eps/mu))^2] (K H) - sinh(K H) = 0
  =>  Etilde_c sqrt(eps/mu) = sqrt( 2( sinh(KH)/(KH) - 1 ) ),

with K = m*pi/l1 the discrete in-plane wavenumber and H = l2 the thickness, so
K H = m*pi/(l1/l2).  The slender limit recovers Euler's column,
Etilde_c sqrt(eps/mu) -> (m*pi)/( (l1/l2) sqrt(3) ).

Writes paper_draft/pics/lam1_elec.pdf .
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import plotstyle as ps
ps.apply_style()
import matplotlib.pyplot as plt

PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")


def Ec_exact(AR, m):
    """Eq. 72: normalized critical nominal field vs aspect ratio AR=l1/l2."""
    KH = m * np.pi / AR
    val = 2.0 * (np.sinh(KH) / KH - 1.0)
    return np.sqrt(val)


def Ec_euler(AR, m):
    """slender (Euler) limit."""
    return m * np.pi / (AR * np.sqrt(3.0))


def main():
    AR = np.linspace(0.6, 10.0, 400)          # aspect ratio l1/l2, high resolution
    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    for i, m in enumerate((1, 2, 3)):
        ax.plot(AR, Ec_exact(AR, m), color=ps.PALETTE[i], label=fr"$m={m}$ (exact)")
        ax.plot(AR, Ec_euler(AR, m), color=ps.PALETTE[i], ls=(0, (5, 4)), lw=0.9)
    ax.set_xlabel(r"aspect ratio $L/H_f$")
    ax.set_ylabel(r"critical nominal field $\tilde E_c\,\sqrt{\epsilon/\mu}$")
    ax.set_xlim(0.6, 10.0)
    ax.set_ylim(0, 6)
    ax.set_title(r"Free-standing block, $\lambda_{pre}=1$ (after Yang \emph{et al.})")
    # legend: solid = exact (Eq.72), dashed = Euler
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=ps.PALETTE[i], label=fr"$m={m}$")
               for i, m in enumerate((1, 2, 3))]
    handles.append(Line2D([0], [0], color="0.35", ls=(0, (5, 4)), lw=0.9,
                          label=r"Euler limit"))
    ax.legend(handles=handles, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(f"{PICS}/lam1_elec.pdf")
    plt.close(fig)
    for m in (1, 2, 3):
        print(f"m={m}: Ec(AR=1)={Ec_exact(1.0,m):.3f}, Ec(AR=5)={Ec_exact(5.0,m):.3f}")
    print("lam1_elec done")


if __name__ == "__main__":
    main()
