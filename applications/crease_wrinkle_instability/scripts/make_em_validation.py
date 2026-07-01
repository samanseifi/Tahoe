#!/usr/bin/env python3
"""Pre-strain dependence at gbar=2: linear theory (line) vs refined-mesh FE
onsets (markers), swept lambda=0.8..1.3.  Reads fea_em/lam_onset.json.
(a) nominal critical field Etilde_c vs pre-strain.
(b) true field E2c = lambda*Etilde_c -- nearly pre-strain-independent."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import plotstyle as ps
ps.apply_style()
import matplotlib.pyplot as plt
import elec_stability as es

ROOT = os.path.join(os.path.dirname(__file__), "..")
PICS = os.path.join(ROOT, "paper_draft", "pics")
FEJSON = os.path.join(ROOT, "fe_tests", "common", "lam_onset.json")
H = 4.0
GBAR = 2.0

# FE (crease) onsets
fe_lam, fe_Ec = [], []
if os.path.exists(FEJSON):
    for r in json.load(open(FEJSON)):
        if r["Ec"] is not None:
            fe_lam.append(r["lam"]); fe_Ec.append(r["Ec"])
fe_lam, fe_Ec = np.array(fe_lam), np.array(fe_Ec)
fe_eps = 1.0 - fe_lam

# theory (wrinkle) at gbar=2
Ks = np.linspace(0.05, 2.6, 80)
lam_th = np.linspace(0.78, 1.32, 18)
th_Ec = np.array([(es.critical_field(l, GBAR, H, Ks) or [np.nan])[0] for l in lam_th])
eps_th = 1.0 - lam_th

fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.7))
# (a) nominal field
ax[0].plot(eps_th, th_Ec, "-", color=ps.PALETTE[0], label="linear theory (wrinkle)")
if len(fe_lam):
    ps.fe_marker(ax[0], fe_eps, fe_Ec, color=ps.PAIR[1], marker="o", label="FE (crease)")
ax[0].set_xlabel(r"pre-strain $\varepsilon^{pre}=1-\lambda_{pre}$")
ax[0].set_ylabel(r"nominal field $\tilde E_c\,\sqrt{\epsilon/\mu}$")
ax[0].set_title(r"Nominal critical field, $\bar\gamma=2$")
ax[0].legend(frameon=False)
# (b) true field
ax[1].plot(eps_th, th_Ec * lam_th, "-", color=ps.PALETTE[0], label="linear theory (wrinkle)")
if len(fe_lam):
    ps.fe_marker(ax[1], fe_eps, fe_Ec * fe_lam, color=ps.PAIR[1], marker="o", label="FE (crease)")
ax[1].set_xlabel(r"pre-strain $\varepsilon^{pre}=1-\lambda_{pre}$")
ax[1].set_ylabel(r"true field $E_{2,c}\,\sqrt{\epsilon/\mu}=\lambda_{pre}\tilde E_c\,\sqrt{\epsilon/\mu}$")
ax[1].set_title(r"True critical field (nearly $\lambda$-independent)")
ax[1].legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{PICS}/em_validation.pdf")
plt.close(fig)
if len(fe_lam):
    print("FE true field E2c=lam*Ec:", np.round(fe_Ec * fe_lam, 3))
print("em_validation.pdf done")
