#!/usr/bin/env python3
"""Tensile-test post-processor — verify J2 cubic-spline hardening (#49).

Reads tensile.io0.exo and extracts the (α, σ_zz, VM_Kirch) history at a
representative IP in the centre of the bar (axis x=y=0, mid-length).  At
that point the bar is in clean uniaxial tension so:

  - σ_VM == |σ_zz - σ_yy| ≈ σ_zz (since σ_xx = σ_yy ≈ 0 free-boundary)
  - The current yield surface is σ_VM = σ_y(α)
  - Plotting σ_VM vs α should overlay the input cubic-spline knots
    to within the radial-return tolerance.

Outputs:
  tensile_results.csv          — per-frame (frame, t, α, σ_zz, VM, p)
  tensile_hardening.png        — VM_Kirch vs α with the input knots overlaid
  tensile_pdelta.png           — total axial reaction vs prescribed δ
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# Spline knots — must match tensile.xml exactly.
SPLINE = [
    (0.000, 250.0), (0.005, 320.0), (0.020, 420.0),
    (0.050, 500.0), (0.080, 535.0), (0.110, 560.0),
    (0.150, 580.0), (0.200, 600.0), (0.300, 625.0),
    (0.400, 640.0), (0.500, 650.0),
]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    stem = sys.argv[1] if len(sys.argv) > 1 else "tensile"
    io0  = os.path.join(here, f"{stem}.io0.exo")
    if not os.path.exists(io0):
        print(f"missing: {io0}"); sys.exit(1)

    f = Dataset(io0, "r")
    nv = [b''.join(r).decode().rstrip('\x00').strip()
          for r in f.variables["name_nod_var"][:]]
    var = {n: f.variables[f"vals_nod_var{i+1}"] for i, n in enumerate(nv)}

    cx, cy, cz = (f.variables[c][:] for c in ("coordx", "coordy", "coordz"))
    nframes = f.dimensions["time_step"].size
    t_whole = f.variables["time_whole"][:]

    # Probe IP at the centre of the bar: choose the node closest to
    # (x=0, y=0, z=Lz/2) where Lz = 5 mm here.
    probe = int(np.argmin(np.abs(cz - 2.5) + np.abs(cx) + np.abs(cy)))
    print(f"probe node: idx={probe}  at (x={cx[probe]:.3f}, y={cy[probe]:.3f}, z={cz[probe]:.3f})")

    rows = []
    for k in range(nframes):
        alpha = float(var["alpha"][k][probe])
        s33   = float(var["s33"][k][probe])
        s11   = float(var["s11"][k][probe])
        s22   = float(var["s22"][k][probe])
        vmk   = float(var["VM_Kirch"][k][probe])
        prs   = float(var["press"][k][probe])
        rows.append((k, float(t_whole[k]), alpha, s33, s11, s22, vmk, prs))

    csv = os.path.join(here, f"{stem}_results.csv")
    with open(csv, "w") as fh:
        fh.write("frame, t, alpha, sigma_zz_MPa, sigma_xx, sigma_yy, VM_Kirch, press\n")
        for r in rows:
            fh.write(", ".join(f"{x:.5f}" for x in r) + "\n")
    print(f"wrote {csv}")

    # Total axial reaction = sum of s_zz weighted by area at the top face.
    # Easier proxy: |s33|·area_xy0 at z=Lz_top.  Not needed for the main
    # validation; just useful as a sanity P-δ.
    top_mask = np.abs(cz - cz.max()) < 1e-9
    Lx_half = cx[top_mask].max()
    Ly_half = cy[top_mask].max()
    area_top = Lx_half * Ly_half     # quarter-symmetry → quarter area
    Ps   = np.array([np.mean(var["s33"][k][top_mask]) * area_top for k in range(nframes)])
    Dzs  = np.array([np.mean(var["D_Z"][k][top_mask]) for k in range(nframes)])

    # --- Hardening overlay ---
    alphas_sim = np.array([r[2] for r in rows])
    VMs        = np.array([r[6] for r in rows])
    s33s       = np.array([r[3] for r in rows])

    a_spline = np.array([p[0] for p in SPLINE])
    s_spline = np.array([p[1] for p in SPLINE])

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(a_spline, s_spline, "o-", color="tab:blue",
            label="input spline knots (XML)", lw=1.5, ms=8)
    ax.plot(alphas_sim, VMs, "s", color="tab:red", ms=6,
            label="Tahoe simulation (VM_Kirch at centre)")
    ax.plot(alphas_sim, s33s, "x", color="tab:green", ms=6,
            label="Tahoe simulation (σ_zz at centre)")
    ax.set_xlabel("equivalent plastic strain  α")
    ax.set_ylabel("yield stress  σ_y  [MPa]")
    ax.set_title(f"J2 hardening verification — {stem}")
    ax.grid(alpha=0.3); ax.legend()
    pfile = os.path.join(here, f"{stem}_hardening.png")
    plt.savefig(pfile, dpi=120, bbox_inches="tight")
    print(f"wrote {pfile}")

    # --- P-δ ---
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(Dzs, Ps, "o-", color="tab:red")
    ax.set_xlabel("top-face displacement  δ_z [mm]")
    ax.set_ylabel("mean axial stress × area  ≈  P_quarter [N]")
    ax.set_title(f"P-δ — {stem}")
    ax.grid(alpha=0.3)
    pfile = os.path.join(here, f"{stem}_pdelta.png")
    plt.savefig(pfile, dpi=120, bbox_inches="tight")
    print(f"wrote {pfile}")

    # Console summary — error vs nearest spline knot at last frame
    last = rows[-1]
    a_last, vmk_last = last[2], last[6]
    nearest_knot = SPLINE[int(np.argmin(np.abs(a_spline - a_last)))]
    err_pct = 100.0 * abs(vmk_last - nearest_knot[1]) / nearest_knot[1]
    print(f"\nSummary (last frame, t={last[1]:.2f}):")
    print(f"  α      = {a_last:.4f}")
    print(f"  σ_zz   = {last[3]:.1f} MPa")
    print(f"  VM     = {vmk_last:.1f} MPa")
    print(f"  nearest knot: ({nearest_knot[0]:.3f}, {nearest_knot[1]:.1f}) MPa")
    print(f"  |VM − σ_y(α)| / σ_y(α) ≈ {err_pct:.2f} %")


if __name__ == "__main__":
    main()
