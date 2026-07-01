#!/usr/bin/env python3
"""No-pre-strain (lambda=1) critical field vs elastocapillary number, gbar up to 20:
linear perturbation theory (smooth wrinkle) against the nonlinear FE onsets
(subcritical crease) from the refined sweep (fea_em/gamma_onset.json)."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import plotstyle as ps
ps.apply_style()
import matplotlib.pyplot as plt
import elec_stability as es

PICS = os.path.join(os.path.dirname(__file__), "..", "paper_draft", "pics")
FEJSON = os.path.join(os.path.dirname(__file__), "..", "fe_tests", "2_em_no_prestrain", "gamma_onset.json")

# theory (wrinkle), lambda=1, gbar 0..20.  At lambda_pre=1 the nominal and true
# fields coincide; we report the NOMINAL field Etilde_c (= Phi_c/H_f), which is
# what the FE applies and measures.
Ks = np.linspace(0.02, 3.0, 1500)
gb = np.linspace(0.0, 20.0, 80)
th = np.array([(es.critical_field(1.0, g, 4.0, Ks) or [np.nan])[0] for g in gb])

# FE (crease) onsets from the refined sweep
fe_gb, fe_Ec = [], []
if os.path.exists(FEJSON):
    for r in json.load(open(FEJSON)):
        if r["Ec"] is not None:
            fe_gb.append(r["gbar"]); fe_Ec.append(r["Ec"])
fe_gb, fe_Ec = np.array(fe_gb), np.array(fe_Ec)

fig, ax = plt.subplots(figsize=(5.4, 4.0))
ax.plot(gb, th, "-", color=ps.PALETTE[1], label=r"linear theory (wrinkle)")
if len(fe_gb):
    ps.fe_marker(ax, fe_gb, fe_Ec, color=ps.PAIR[1], marker="o",
                 label=r"finite element (crease)")
ax.set_xlabel(r"elastocapillary number $\gamma/(\mu H_f)$")
ax.set_ylabel(r"nominal critical field $\tilde E_c\,\sqrt{\epsilon/\mu}$ at $\lambda_{pre}=1$")
ax.set_xlim(0, 20); ax.set_ylim(0.0, 5.0)
ax.legend(frameon=False, fontsize=9)
ax.set_title(r"No pre-strain: smooth wrinkle (theory) vs subcritical crease (FE)")
fig.tight_layout()
fig.savefig(f"{PICS}/lam1_compare.pdf")
plt.close(fig)
if len(fe_gb):
    print("FE/theory ratio:", np.round(fe_Ec/np.interp(fe_gb, gb, th), 3))
print("lam1_compare.pdf done")
