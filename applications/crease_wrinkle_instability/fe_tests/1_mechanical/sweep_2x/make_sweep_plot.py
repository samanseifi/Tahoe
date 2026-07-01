#!/usr/bin/env python3
"""Plot elastocapillary sweep: critical strain eps_c and wrinkle wavelength l/H vs
gbar -- FE (sweep_results.json) against the accurate linear-perturbation theory
(nlayer_stability.critical, the same solver behind Fig 4 gl_mech/ge_mech).
 - eps_c panel extends to gbar=5 and meets the Biot dashed line as gbar->0.
 - wavelength panel shows BOTH FE estimates (FFT peak, autocorrelation); gbar=0 is
   omitted (scale-free Biot problem -> no selected wavelength)."""
import os, sys, json
import numpy as np
SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "..", "scripts")
sys.path.insert(0, SCRIPTS)
import plotstyle as ps
ps.apply_style()
import matplotlib.pyplot as plt
import nlayer_stability as nl

HERE = os.path.dirname(os.path.abspath(__file__))
H = 8.0
BIOT = 0.456

# FE results
fe = json.load(open(os.path.join(HERE, "sweep_results.json")))
fg = np.array([r["gbar"] for r in fe])
fe_eps = np.array([r["eps_c"] if r["eps_c"] is not None else np.nan for r in fe])
fe_lf = np.array([r["lH_fft"] if r.get("lH_fft") is not None else np.nan for r in fe])
fe_la = np.array([r["lH_ac"] if r.get("lH_ac") is not None else np.nan for r in fe])

# accurate theory (nlayer_stability, 300 log-spaced wavenumbers)
Ks = np.geomspace(0.04, 6.0, 200)
gb_e = np.linspace(0.0, 5.0, 36)
th_e = np.array([nl.critical([1.0], [H], [], gb * H, Ks, base="sliding")[0] for gb in gb_e])
gb_l = np.linspace(0.1, 5.0, 36)                 # skip gbar=0 (scale-free)
th_l = np.array([2 * np.pi / nl.critical([1.0], [H], [], gb * H, Ks, base="sliding")[1] / H for gb in gb_l])

fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8))

# (a) critical strain
ax[0].plot(gb_e, th_e, "-", color=ps.PALETTE[1], label="linear theory")
ps.fe_marker(ax[0], fg, fe_eps, color=ps.PAIR[1], marker="o", label="FE (symmetric)")
ax[0].axhline(BIOT, color="0.6", lw=0.9, ls="--")
ax[0].text(3.4, BIOT + 0.008, r"Biot ($\bar\gamma\!\to\!0$)", fontsize=8, color="0.4")
ax[0].set_xlabel(r"elastocapillary number $\gamma/(\mu H_f)$")
ax[0].set_ylabel(r"critical strain $\varepsilon_c$")
ax[0].set_xlim(0, 5); ax[0].set_title("Onset of instability")
ax[0].legend(frameon=False)

# (b) wavelength -- two FE methods vs theory
ax[1].plot(gb_l, th_l, "-", color=ps.PALETTE[1], label="linear theory")
ps.fe_marker(ax[1], fg, fe_lf, color=ps.PAIR[1], marker="o", label="FE (FFT)")
ps.fe_marker(ax[1], fg, fe_la, color=ps.PALETTE[3], marker="s", label="FE (autocorr)")
ax[1].set_xlabel(r"elastocapillary number $\gamma/(\mu H_f)$")
ax[1].set_ylabel(r"wavelength $\ell_c/H_f$")
ax[1].set_xlim(0, 5); ax[1].set_title("Wrinkle wavelength")
ax[1].legend(frameon=False)

fig.tight_layout(); fig.savefig(os.path.join(HERE, "sweep_eps_wavelength.pdf"))
plt.close(fig)
print("sweep_eps_wavelength.pdf written")
for r in fe:
    print(f"  gbar={r['gbar']}: eps_c={r['eps_c']} lH_fft={r.get('lH_fft')} lH_ac={r.get('lH_ac')}")
