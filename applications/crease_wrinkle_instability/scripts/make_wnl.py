#!/usr/bin/env python3
"""Weakly nonlinear figure: Landau coefficient LAM2 vs elastocapillary number,
from direct nonlinear continuation (nl_continuation.py). eps - eps_c = -LAM2 A^2.
LAM2>0 subcritical (crease-like); -> 0 near gbar~2 (subcriticality killed by ST)."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import plotstyle as ps
ps.apply_style()
import matplotlib.pyplot as plt

PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")
d = json.load(open(os.path.join(os.path.dirname(__file__), "..", "fe_tests", "common", "wnl_landau.json")))
gb = np.array([r["gbar"] for r in d]); L2 = np.array([r["LAM2"] for r in d])

fig, ax = plt.subplots(figsize=(5.4, 4.0))
ax.axhline(0.0, color="0.6", lw=0.8)
ax.plot(gb, L2, "-", color=ps.PALETTE[1])
ps.fe_marker(ax, gb, L2, color=ps.PAIR[1], marker="o", label=r"continuation")
ax.fill_between([0, 2.0], [-0.02, -0.02], [0.13, 0.13], color=ps.PALETTE[2], alpha=0.06)
ax.text(0.6, 0.085, "subcritical\n(crease-like)", color=ps.PALETTE[2], fontsize=9, ha="center")
ax.text(2.5, 0.03, "marginal /\nsupercritical", color=ps.PALETTE[1], fontsize=9, ha="center")
ax.set_xlabel(r"elastocapillary number $\gamma/(\mu H_f)$")
ax.set_ylabel(r"Landau coefficient $\Lambda_2$  ($\varepsilon-\varepsilon_c=-\Lambda_2 A^2$)")
ax.set_xlim(0, 3.1); ax.set_ylim(-0.02, 0.13)
ax.set_title(r"Surface tension weakens the subcritical wrinkle")
fig.tight_layout(); fig.savefig(f"{PICS}/wnl_landau.pdf"); plt.close(fig)
print("wnl_landau.pdf done")
