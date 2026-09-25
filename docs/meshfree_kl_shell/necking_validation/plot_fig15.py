#!/usr/bin/env python3
"""Overlay a necking run on the digitized Wang & Bazilevs Fig 15: effective stress |R_z|/A0 vs U_norm.
Usage: plot_fig15.py run.log [label]   (A0 = 2 pi R h with R=10, h=1)"""
import csv, math, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fig18_validation", "reference_data", "fig15_necking_effstress_vs_Unorm.csv")
A0 = 2*math.pi*10.0*1.0
log = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else "Tahoe RKShellT"
P = [l.split() for l in open(log) if l.startswith("[RKShell-fig15]")]
U = [float(p[1]) for p in P]; S = [abs(float(p[4]))/A0 for p in P]
ref = list(csv.DictReader(open(REF)))
fig, ax = plt.subplots(figsize=(7, 4.5)); INK2 = "#52514e"
for key, st, lab in [("rkpm_membrane_3pt", "-", "paper RKPM 3PT + membrane stab"), ("nostab_3pt", "--", "paper 3PT no stab"),
                     ("ambati_2018", ":", "Ambati 2018"), ("alaydin_2021", "-.", "Alaydin 2021")]:
    pts = [(float(r["U_norm_mm"]), float(r[key])) for r in ref if r[key]]
    ax.plot(*zip(*pts), st, color=INK2, lw=1.4, label=lab)
ax.plot(U, S, color="#1baf7a", lw=2, label=label)
ax.set_xlabel("U_norm [mm]"); ax.set_ylabel("effective stress |R|/A0 [MPa]"); ax.set_xlim(0, 12.5); ax.set_ylim(0, 700)
ax.grid(True, color="#e4e3df"); ax.legend(frameon=False, fontsize=8, loc="lower right"); fig.tight_layout(); fig.savefig("fig15_compare.png", dpi=150)
peak = max(S) if S else 0
print("points %d  last U_norm %.3f  peak stress %.1f MPa (paper peak ~609)" % (len(U), U[-1] if U else 0, peak))
for u0 in (0.5, 1, 2, 4, 6, 8, 10, 12):
    r = [rr for rr in ref if abs(float(rr["U_norm_mm"]) - u0) < 0.13]
    if not r or not U: continue
    i = min(range(len(U)), key=lambda k: abs(U[k]-u0)); mine = S[i] if abs(U[i]-u0) < 0.3 else float("nan")
    print("U=%5.1f  paper %6s  ambati %6s  ours %6.1f" % (u0, r[0]["rkpm_membrane_3pt"], r[0]["ambati_2018"], mine))
