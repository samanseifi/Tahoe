#!/usr/bin/env python3
"""Render the deformed meshfree cube (90-degree twist) from the ExodusII output."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from netCDF4 import Dataset

d = Dataset("cube_rotate.io0.exo")
x = d.variables['coordx'][:]; y = d.variables['coordy'][:]; z = d.variables['coordz'][:]
names = [bytes(c).split(b'\x00')[0].decode() for c in d.variables['name_nod_var'][:]]
def var(nm, t): i = names.index(nm)+1; return d.variables[f'vals_nod_var{i}'][t,:]
nt = d.variables['time_whole'].shape[0]
times = d.variables['time_whole'][:]

fig = plt.figure(figsize=(13, 4.2))
frames = [0, nt//2, nt-1]
for p, t in enumerate(frames):
    ux, uy, uz = var('D_X', t), var('D_Y', t), var('D_Z', t)
    X, Y, Z = x+ux, y+uy, z+uz
    twist = np.degrees(np.arctan2(Z, Y) - np.arctan2(z, y))
    twist = (twist + 360) % 360
    twist[np.hypot(y, z) < 1e-9] = 0
    ax = fig.add_subplot(1, 3, p+1, projection='3d')
    s = ax.scatter(X, Y, Z, c=twist, cmap='viridis', s=18, vmin=0, vmax=90)
    ax.set_title(f"t = {times[t]:.2f}   (twist {times[t]*90:.0f}°)")
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_xlim(0, 1); ax.set_ylim(-0.7, 0.7); ax.set_zlim(-0.7, 0.7)
    ax.view_init(elev=18, azim=-65)
fig.colorbar(s, ax=fig.axes, shrink=0.6, label='twist angle (deg)', pad=0.02)
fig.suptitle("Meshfree (RKPM) cube, Neo-Hookean — fixed face + 90° twist", fontsize=12)
fig.savefig("cube_twist.png", dpi=130, bbox_inches='tight')
print("wrote cube_twist.png")

# twist-vs-axis profile at final step
ux, uy, uz = var('D_X', -1), var('D_Y', -1), var('D_Z', -1)
X, Y, Z = x+ux, y+uy, z+uz
col = [i for i in range(len(x)) if abs(y[i]-0.5) < 1e-6 and abs(z[i]-0.5) < 1e-6]
col.sort(key=lambda i: x[i])
xs = [x[i] for i in col]
tw = [((np.degrees(np.arctan2(Z[i], Y[i]) - np.arctan2(z[i], y[i])) + 360) % 360) for i in col]
fig2, ax2 = plt.subplots(figsize=(5, 3.6))
ax2.plot(xs, tw, 'o-', label='meshfree result')
ax2.plot([0, 1], [0, 90], 'k--', lw=1, label='linear (rigid torsion)')
ax2.set_xlabel('axial position x'); ax2.set_ylabel('twist angle (deg)')
ax2.set_title('Twist along the axis (corner fiber)'); ax2.legend(); ax2.grid(alpha=0.3)
fig2.savefig("cube_twist_profile.png", dpi=130, bbox_inches='tight')
print("wrote cube_twist_profile.png")
