#!/usr/bin/env python3
"""Plot the elasto-plastic pinched-cylinder load-displacement curve (Fig 18 analog)."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = np.loadtxt("fig18_curve.txt")
if d.ndim == 1: d = d.reshape(1, -1)
u = -d[:,0]            # crush displacement (mm, positive inward)
R = -d[:,1]            # reaction (positive = resisting)
fig, ax = plt.subplots(figsize=(6,4.2))
ax.plot(u, R, '-o', ms=3, lw=1.4, color="#1f4e9c")
ax.set_xlabel("load-point crush displacement  |u_x|  (mm)")
ax.set_ylabel("reaction force  |R_x|")
ax.set_title("Elasto-plastic pinched cylinder crush (meshfree KL shell)")
ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig("fig18.png", dpi=130)
print(f"wrote fig18.png  ({len(u)} points, u_max={u.max():.1f}mm, R_max={R.max():.3e})")
