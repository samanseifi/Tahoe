#!/usr/bin/env python3
"""Plot the elasto-plastic pinched-cylinder load-displacement curve (Fig 18 analog).

Reads fig18_hardening_curve.txt (the corrected nt=40 quasi-static crush). Shows the raw
reaction (explicit dynamics retains some kinetic-energy oscillation) plus a moving-average
trend revealing the underlying elastic -> yield -> plastic response.
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(fn):
    if not os.path.exists(fn): return None
    d = np.loadtxt(fn)
    if d.size == 0: return None
    if d.ndim == 1: d = d.reshape(1, -1)
    return -d[:,0], -d[:,1]

def smooth(y, w=5):
    if len(y) < w: return y
    yp = np.pad(y, w//2, mode="edge")          # edge-pad so endpoints aren't pulled down
    return np.convolve(yp, np.ones(w)/w, mode="valid")[:len(y)]

fig, ax = plt.subplots(figsize=(6.6,4.5))
r = load("fig18_hardening_curve.txt") or load("fig18_curve.txt")
if r is not None:
    u, R = r
    ax.plot(u, R, '-', lw=0.8, color="#9bb8e0", alpha=0.9, label="reaction (raw, explicit)")
    ax.plot(u, smooth(R), '-', lw=2.2, color="#1f4e9c", label="moving-average trend")
ax.set_xlabel("load-point crush displacement  |u_x|  (mm)")
ax.set_ylabel("reaction force  |R_x|")
ax.set_title("Elasto-plastic pinched cylinder crush (meshfree KL shell, nt=160 (8320 nodes))")
ax.grid(True, alpha=0.3); ax.legend(loc="lower right", fontsize=9)
fig.tight_layout(); fig.savefig("fig18.png", dpi=130)
print("wrote fig18.png")
