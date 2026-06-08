#!/usr/bin/env python3
"""Plot the elasto-plastic pinched-cylinder load-displacement curve(s) (Fig 18 analog).

Reads fig18_curve.txt (always) and, if present, fig18_hardening_curve.txt /
fig18_perfect.txt to overlay the hardening vs near-perfect-plasticity responses.
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(fn):
    if not os.path.exists(fn): return None
    d = np.loadtxt(fn)
    if d.size == 0: return None
    if d.ndim == 1: d = d.reshape(1, -1)
    return -d[:,0], -d[:,1]      # crush displacement (mm), resisting reaction

fig, ax = plt.subplots(figsize=(6.4,4.4))
series = [("fig18_hardening_curve.txt", "linear hardening  (H=3e5)", "#1f4e9c", "-o"),
          ("fig18_perfect.txt",         "near-perfect plasticity  (H=0)", "#c0392b", "-s")]
plotted = False
for fn, lab, col, st in series:
    r = load(fn)
    if r is None: continue
    ax.plot(r[0], r[1], st, ms=3, lw=1.5, color=col, label=lab); plotted = True
if not plotted:
    r = load("fig18_curve.txt")
    if r is not None: ax.plot(r[0], r[1], '-o', ms=3, color="#1f4e9c", label="crush"); plotted = True

ax.set_xlabel("load-point crush displacement  |u_x|  (mm)")
ax.set_ylabel("reaction force  |R_x|")
ax.set_title("Elasto-plastic pinched cylinder crush (meshfree KL shell)")
ax.grid(True, alpha=0.3); ax.legend(loc="lower right", fontsize=9)
fig.tight_layout(); fig.savefig("fig18.png", dpi=130)
print("wrote fig18.png")
