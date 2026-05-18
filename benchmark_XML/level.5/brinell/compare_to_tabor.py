#!/usr/bin/env python3
"""Brinell post-processor — extract P-δ and the Tabor relation HB/σ_y (#47).

Reads <stem>.io0.exo (block 1 — indenter) and <stem>.io1.exo (block 2 —
the plastic block / contact group), pairs each time frame with the
prescribed top-face displacement δ from the time history, and computes:

  - δ(t)        — indentation depth (apex displacement of indenter)
  - P(t)        — total contact load on the quarter (full sphere = 4P)
  - a(t)        — contact radius (largest r on a striker with F_z > 5 % of peak)
  - p_m(t)      — mean contact pressure  P / (π a²)
  - p_m / σ_y0  — Tabor ratio (should → ~2.8 in the fully-plastic regime)
  - HB          — Brinell hardness, computed as P / (π D h) where D = 2R
                   and h is the impression depth from the spherical cap

Writes:
  brinell_results.csv  — one row per output frame
  brinell_pdelta.png   — P-δ curve with elastic Hertz reference overlaid
  brinell_tabor.png    — p_m / σ_y0 vs δ/R, with Tabor line at 2.8

Usage:
  python3 compare_to_tabor.py [brinell|brinell_smoke]   # default: brinell_smoke
"""

import os
import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# --- Constants matching brinell{,_smoke}.xml ----------------------------------
R       = 5.0          # indenter radius, mm
E_BLOCK = 200_000.0    # block Young's modulus, MPa
NU      = 0.3
SIG_Y0  = 250.0        # block initial yield, MPa
# Effective contact modulus for one rigid + one elastic body:
#   1/E* = (1-ν²)/E_block   (rigid → infinite E)
ESTAR_RIGID = E_BLOCK / (1.0 - NU * NU)


def striker_areas_from_geom(geom_path):
    """Tributary striker areas from the side-set 1 (indenter bottom, ns_id=7).

    Same algorithm as level.5/hertz/compare_to_analytical.py: each quad
    facet contributes a quarter of its area to each of its four nodes.
    """
    with open(geom_path) as f:
        text = f.read()

    nodes_idx    = text.rfind("\n*nodes\n") + 1
    sidesets_idx = text.index("*sidesets")
    elements_idx = text.index("*elements")

    lines = text[nodes_idx:].split("\n")
    n_nod = int(lines[1].strip())
    coords = {}
    for i in range(3, 3 + n_nod):
        parts = lines[i].split()
        coords[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))

    blocks = text[elements_idx:nodes_idx].split("*set")[1:]
    elem_conn = {}
    for bi, blk in enumerate(blocks):
        ll = blk.strip().split("\n")
        nelem = int(ll[0].split()[0])
        nen   = int(ll[1].split()[0])
        for k in range(nelem):
            parts = ll[2 + k].split()
            elem_conn[(bi + 1, int(parts[0]))] = list(map(int, parts[1:1 + nen]))

    ss = text[sidesets_idx:elements_idx].split("*set")[1:]
    ss1 = ss[0].strip().split("\n")
    n_ss1 = int(ss1[0])
    pairs = [tuple(map(int, ln.split())) for ln in ss1[1:1 + n_ss1]]

    # Hex8 face → local node ordering (matches Tahoe HexT convention)
    face_local = {1: [0, 3, 2, 1], 2: [4, 5, 6, 7]}
    striker_area = {}
    for eid, fc in pairs:
        ns = elem_conn[(1, eid)]
        nf = [ns[k] for k in face_local[fc]]
        p = [np.array(coords[n]) for n in nf]
        face_area = (0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
                     + 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[0])))
        for n in nf:
            striker_area[n] = striker_area.get(n, 0.0) + 0.25 * face_area
    return striker_area


def analyse_frame(f0, f1, frame):
    """Single time-step extraction."""
    # io0 contains all nodes from BOTH blocks (indenter + plastic block).
    # var1..3 = D_X, D_Y, D_Z;  var4..9 = s11..s12 (Cauchy stress).
    DZ0 = f0.variables["vals_nod_var3"][frame]
    cz0 = f0.variables["coordz"][:]
    cx0 = f0.variables["coordx"][:]
    cy0 = f0.variables["coordy"][:]

    # δ = indentation depth.  For a rigid indenter on a half-space the
    # cleanest unambiguous definition is the downward displacement of the
    # indenter top face — equal to the prescribed BC at NS1.  In our quarter
    # mesh the indenter top apex is the node with the largest cz on the
    # x = y = 0 column.
    on_axis = (np.abs(cx0) < 1e-9) & (np.abs(cy0) < 1e-9)
    cz_axis = np.where(on_axis, cz0, -1.0e18)
    top_apex = int(np.argmax(cz_axis))
    delta = -float(DZ0[top_apex])      # downward → positive δ

    # io1 — contact pair output.  vals_nod_var1..3 are D_X,Y,Z;
    # vals_nod_var6 is F_z on the striker (when stress=1 isn't requested,
    # only displacements are written.  See solid_element_nodal_output in XML).
    # For PenaltyContact3DT, only displacements are output at strikers unless
    # we add stress=1 there too — this is a known limitation; we fall back
    # to integrating the indenter bottom-face reaction force via
    # ∫_S σ_n dA from the block's surface stress instead.
    #
    # Simpler, robust route: read the block (io1.exo) and integrate F_z
    # at the contact patch surface using the striker-area weights against
    # the *equivalent reaction force* derived from displacement compatibility.
    # Here we use the block's vals_nod_var with whatever is available —
    # if F_z isn't written, return NaN and the user can re-run with
    # appropriate output flags.
    try:
        Fz = f1.variables["vals_nod_var6"][frame]
    except (KeyError, IndexError):
        Fz = None

    if Fz is None:
        # No striker-force output: compute P from displacement BC reaction at
        # node-set 1 (indenter top).  This needs a separate output channel we
        # haven't requested; for now return NaN-padded result.
        return {"delta": delta, "P": np.nan, "a": np.nan, "p_m": np.nan}

    nmap = f1.variables["node_num_map"][:]
    cx1, cy1 = cx0[nmap - 1], cy0[nmap - 1]
    DX1 = f0.variables["vals_nod_var1"][frame][nmap - 1]
    DY1 = f0.variables["vals_nod_var2"][frame][nmap - 1]
    r = np.sqrt((cx1 + DX1)**2 + (cy1 + DY1)**2)
    P = float(np.sum(Fz))                          # quarter
    in_contact = Fz > 0.05 * max(Fz.max(), 1e-12)
    a = float(r[in_contact].max()) if in_contact.any() else 0.0
    p_m = P / (math.pi * a * a) if a > 0 else 0.0
    return {"delta": delta, "P": P, "a": a, "p_m": p_m}


def hertz_P(delta, R, Estar):
    """Hertz quarter-sphere load: P = (4/3) E* sqrt(R) δ^{3/2} / 4."""
    if delta <= 0: return 0.0
    return (4.0 / 3.0) * Estar * math.sqrt(R) * delta**1.5 / 4.0


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    stem = sys.argv[1] if len(sys.argv) > 1 else "brinell_smoke"
    io0  = os.path.join(here, f"{stem}.io0.exo")   # indenter (group 1)
    io1  = os.path.join(here, f"{stem}.io1.exo")   # plastic block (group 2, has alpha/VM_Kirch/press)
    io2  = os.path.join(here, f"{stem}.io2.exo")   # contact group (F_X/Y/Z on strikers)
    geom = os.path.join(here, f"{stem}.geom")
    for p in (io0, io2, geom):
        if not os.path.exists(p):
            print(f"missing: {p}")
            sys.exit(1)

    f0 = Dataset(io0, "r")
    f2 = Dataset(io2, "r")     # contact group — has F_z on strikers
    n_frames = f0.dimensions["time_step"].size

    rows = []
    for k in range(n_frames):
        r = analyse_frame(f0, f2, k)
        rows.append((k, r["delta"], r["P"], r["a"], r["p_m"]))

    # Write CSV
    csv = os.path.join(here, f"{stem}_results.csv")
    with open(csv, "w") as fh:
        fh.write("frame, delta_mm, P_quarter_N, a_mm, p_mean_MPa, p_m/sigma_y0\n")
        for k, d, P, a, pm in rows:
            ratio = pm / SIG_Y0 if not math.isnan(pm) else float("nan")
            fh.write(f"{k}, {d:.6f}, {P:.4f}, {a:.4f}, {pm:.2f}, {ratio:.3f}\n")
    print(f"wrote {csv}")

    deltas = np.array([r[1] for r in rows])
    Ps     = np.array([r[2] for r in rows])
    p_ms   = np.array([r[4] for r in rows])

    # --- P-δ plot ---
    plt.figure(figsize=(7, 5))
    plt.plot(deltas, Ps, "o-", label="Tahoe (J2 plastic)", color="tab:red")
    d_fine = np.linspace(0, deltas.max(), 50)
    P_hertz = np.array([hertz_P(d, R, ESTAR_RIGID) for d in d_fine])
    plt.plot(d_fine, P_hertz, "--", label="Hertz (elastic-only)", color="tab:blue")
    plt.xlabel("indentation depth δ [mm]")
    plt.ylabel("contact load P (quarter) [N]")
    plt.title(f"Brinell P-δ — {stem}")
    plt.grid(alpha=0.3); plt.legend()
    pfile = os.path.join(here, f"{stem}_pdelta.png")
    plt.savefig(pfile, dpi=120, bbox_inches="tight")
    print(f"wrote {pfile}")

    # --- Tabor plot: p_m / σ_y vs δ/R ---
    plt.figure(figsize=(7, 5))
    plt.plot(deltas / R, p_ms / SIG_Y0, "o-", color="tab:red", label="Tahoe")
    plt.axhline(2.8, ls="--", color="k", alpha=0.5, label="Tabor (2.8 σ_y)")
    plt.xlabel("δ / R")
    plt.ylabel("p_m / σ_y0")
    plt.title(f"Tabor relation — {stem}")
    plt.grid(alpha=0.3); plt.legend()
    tfile = os.path.join(here, f"{stem}_tabor.png")
    plt.savefig(tfile, dpi=120, bbox_inches="tight")
    print(f"wrote {tfile}")

    # --- plastic-strain + von-Mises maps on the block's top surface ---
    # io1 (group 2 — the plastic block) holds 'alpha' (ε_p) and 'VM_Kirch'
    # nodal values.  We slice at z ≈ 0 (the original top face) — those nodes
    # carry the surface impression and yield-plume footprint.
    if os.path.exists(io1):
        f1 = Dataset(io1, "r")
        nv = [b''.join(r).decode().rstrip('\x00').strip()
              for r in f1.variables["name_nod_var"][:]]
        cx, cy, cz = (f1.variables[c][:] for c in ("coordx", "coordy", "coordz"))
        # block-top nodes: z = 0 in the reference config
        top = np.where(np.abs(cz) < 1e-9)[0]
        if "alpha" in nv and "VM_Kirch" in nv and len(top) > 0:
            alpha   = f1.variables[f"vals_nod_var{nv.index('alpha')+1}"][-1][top]
            vmk     = f1.variables[f"vals_nod_var{nv.index('VM_Kirch')+1}"][-1][top]
            DX, DY  = (f1.variables[f"vals_nod_var{nv.index(c)+1}"][-1][top]
                       for c in ("D_X", "D_Y"))
            xT, yT  = cx[top] + DX, cy[top] + DY    # deformed top-face coords

            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            sc0 = axes[0].tricontourf(xT, yT, alpha, levels=20, cmap="plasma")
            axes[0].set_aspect("equal")
            axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
            axes[0].set_title("equivalent plastic strain  α  (block top, last frame)")
            plt.colorbar(sc0, ax=axes[0], label="α")

            sc1 = axes[1].tricontourf(xT, yT, vmk, levels=20, cmap="viridis")
            axes[1].set_aspect("equal")
            axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("y [mm]")
            axes[1].set_title("von Mises Kirchhoff stress [MPa]")
            plt.colorbar(sc1, ax=axes[1], label="σ_VM [MPa]")
            fig.suptitle(f"Block top surface — {stem}")
            ffile = os.path.join(here, f"{stem}_field.png")
            plt.savefig(ffile, dpi=120, bbox_inches="tight")
            print(f"wrote {ffile}")
            print(f"  α range  = [{alpha.min():.4f}, {alpha.max():.4f}]")
            print(f"  σ_VM     = [{vmk.min():.1f}, {vmk.max():.1f}] MPa")

    # Console summary
    print("\nSummary (last frame):")
    last = rows[-1]
    print(f"  δ        = {last[1]:.4f} mm   ({last[1]/R*100:.2f} % of R)")
    print(f"  P        = {last[2]:.2f} N    (quarter; full = {4*last[2]:.1f})")
    print(f"  a        = {last[3]:.4f} mm")
    print(f"  p_m      = {last[4]:.2f} MPa")
    if last[4] > 0:
        print(f"  p_m/σ_y0 = {last[4]/SIG_Y0:.3f}   (Tabor target ≈ 2.8 fully plastic)")


if __name__ == "__main__":
    main()
