#!/usr/bin/env python3
"""
Case 2 figure: critical electric field vs elastocapillary number for a
pre-compressed and a pre-stretched single DE layer, from the corrected
electromechanical bifurcation criterion (elec_stability.det B = 0).
Writes paper_draft/pics/wrinkle_elec1.pdf  (fig:elecfig_plot).

Run only AFTER elec_stability.py --validate confirms the E->0 (Biot) and
Wang-Zhao limits.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
import elec_stability as es

PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")
PAIR = ["#1d3557", "#c75643"]; REF = "#8a8f99"
plt.rcParams.update({"text.usetex": True, "font.family": "serif",
    "text.latex.preamble": r"\usepackage{amsmath}", "font.size": 13,
    "axes.labelsize": 14, "legend.fontsize": 11, "lines.linewidth": 2.0,
    "savefig.bbox": "tight"})


def curve(lam, gbars, H, Ks):
    out = []
    for gb in gbars:
        res = es.critical_field(lam, gb, H, Ks)   # new API takes gbar directly
        out.append(res[0] if res else np.nan)
    return np.array(out)


def main():
    H = 4.0
    Ks = np.linspace(0.05, 4.0, 160)
    gbars = np.linspace(0.0, 5.0, 21)
    Ec_comp = curve(0.8, gbars, H, Ks)   # pre-compression eps_pre=0.2
    Ec_str  = curve(1.2, gbars, H, Ks)   # pre-stretch
    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    ax.plot(gbars, Ec_comp, "-", color=PAIR[0], label=r"pre-compression ($\lambda^{pre}=0.8$)")
    ax.plot(gbars, Ec_str,  "--", color=PAIR[1], label=r"pre-stretch ($\lambda^{pre}=1.2$)")
    ax.set_xlabel(r"elastocapillary number $\gamma/(\mu H_f)$")
    ax.set_ylabel(r"critical field $\tilde E_c\sqrt{\epsilon/\mu}$")
    ax.grid(alpha=0.3); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{PICS}/wrinkle_elec1.pdf"); plt.close(fig)
    print("compression Ec:", np.round(Ec_comp, 2))
    print("stretch     Ec:", np.round(Ec_str, 2))
    print(f"wrote {PICS}/wrinkle_elec1.pdf")


if __name__ == "__main__":
    main()
