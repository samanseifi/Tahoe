#!/usr/bin/env python3
"""Compare Tahoe Hertz contact simulations to the analytical Hertz solution.

Reads the Exodus output of:
  hertz_explicit.io0.exo  — sphere + base nodal field
  hertz_explicit.io1.exo  — contact element strikers (NS7) with F_X,F_Y,F_Z
  hertz_implicit.io0.exo  — implicit run (same convention)
  hertz_implicit.io1.exo

For each:
  1. Identifies the effective indentation δ as the gap closure between
     the sphere apex and the base apex.
  2. Sums quarter-domain reaction force F_quarter (×4 for full).
  3. Determines the contact radius a (largest r where penetration > 0).
  4. Builds p(r) ≈ -F_z(striker) / striker_area, plotted vs analytical
     p(r) = p₀ √(1 - (r/a)²).

Analytical (frictionless, identical materials, half-space):
   E* = E / (2(1-ν²))
   a  = √(R · δ)
   F  = (4/3) E* √R · δ^{3/2}    (full sphere)
   p₀ = 3 F / (2π a²)
"""
import math
import os
import sys

import numpy as np
from netCDF4 import Dataset

R   = 10.0
E   = 200000.0
NU  = 0.3
DELTA_BC = 0.05            # prescribed sphere-top displacement (mm)


def hertz(delta, Estar):
    """Return analytical (a, F_full, p0) given effective indentation delta."""
    a   = math.sqrt(R * delta)
    F   = (4.0 / 3.0) * Estar * math.sqrt(R) * delta ** 1.5
    p0  = 3.0 * F / (2.0 * math.pi * a * a)
    return a, F, p0


def Estar_for(E1, nu1, E2, nu2):
    """Hertz effective contact modulus for two bodies."""
    return 1.0 / ((1.0 - nu1 * nu1) / E1 + (1.0 - nu2 * nu2) / E2)


def striker_areas_from_geom(geom_path):
    """Read the .geom file, extract NS7 (sphere bottom) tributary areas.

    Each striker's area is computed as 1/4 of the sum of areas of the
    quad facets it belongs to (a node touches up to 4 quads in a regular
    grid).  The facets are reconstructed from the SS1 sideset.
    """
    with open(geom_path) as f:
        text = f.read()

    # *nodes is the *last* coordinate section (after *nodesets, *sidesets,
    # *elements).  Find via rfind to skip *nodesets.
    nodes_idx = text.rfind("\n*nodes\n") + 1
    sidesets_idx = text.index("*sidesets")
    elements_idx = text.index("*elements")

    # parse nodes
    nodes_block = text[nodes_idx:]
    lines = nodes_block.split("\n")
    # line 0: "*nodes", line 1: "<n_nod>", line 2: "<ndim>", lines 3..3+n_nod-1: data
    n_nod = int(lines[1].strip())
    coords = {}
    for i in range(3, 3 + n_nod):
        parts = lines[i].split()
        nid = int(parts[0])
        coords[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))

    # parse element blocks (between *elements and *nodes)
    elements_block = text[elements_idx:nodes_idx]
    blocks = elements_block.split("*set")[1:]
    elem_conn = {}  # (block_idx, local_eid) -> [n0..n7]
    for bi, blk in enumerate(blocks):
        ll = blk.strip().split("\n")
        nelem = int(ll[0].split()[0])
        nen = int(ll[1].split()[0])
        for k in range(nelem):
            parts = ll[2 + k].split()
            eid = int(parts[0])
            ns  = list(map(int, parts[1:1 + nen]))
            elem_conn[(bi + 1, eid)] = ns

    # parse SS1 (sphere bottom; block 1, face 1 = nodes {0,3,2,1})
    ss_block = text[sidesets_idx:elements_idx]
    ss = ss_block.split("*set")[1:]
    ss1_lines = ss[0].strip().split("\n")
    ss1_pairs = []
    n_ss1 = int(ss1_lines[0])
    for ln in ss1_lines[1:1 + n_ss1]:
        e, fc = ln.split()
        ss1_pairs.append((int(e), int(fc)))

    # accumulate striker area: each quad face area / 4 to each of its nodes
    face_local = {1: [0, 3, 2, 1], 2: [4, 5, 6, 7]}
    striker_area = {}
    for eid, fc in ss1_pairs:
        ns = elem_conn[(1, eid)]
        face_nodes = [ns[k] for k in face_local[fc]]
        p0 = np.array(coords[face_nodes[0]])
        p1 = np.array(coords[face_nodes[1]])
        p2 = np.array(coords[face_nodes[2]])
        p3 = np.array(coords[face_nodes[3]])
        # quad area = |((p1-p0) x (p2-p0))| / 2 + |((p2-p0) x (p3-p0))| / 2
        a1 = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0))
        a2 = 0.5 * np.linalg.norm(np.cross(p2 - p0, p3 - p0))
        face_area = a1 + a2
        for nid in face_nodes:
            striker_area[nid] = striker_area.get(nid, 0.0) + 0.25 * face_area

    return striker_area, coords


def analyse(io0_path, io1_path, geom_path, label):
    """Read sim outputs and return (delta_eff, F_quarter, a_sim, r_arr, p_arr)."""
    f0 = Dataset(io0_path, "r")
    f1 = Dataset(io1_path, "r")

    cx = f0.variables["coordx"][:]
    cy = f0.variables["coordy"][:]
    cz = f0.variables["coordz"][:]
    DX = f0.variables["vals_nod_var1"][-1]   # last frame
    DY = f0.variables["vals_nod_var2"][-1]
    DZ = f0.variables["vals_nod_var3"][-1]

    # sphere apex: nid 1 (i=j=k=0); base apex: at (0,0,0)
    sphere_apex_idx = 0
    base_apex_idx = np.where(
        (np.abs(cx) < 1e-9) & (np.abs(cy) < 1e-9) & (np.abs(cz) < 1e-9))[0][0]
    sphere_apex_z = cz[sphere_apex_idx] + DZ[sphere_apex_idx]
    base_apex_z  = cz[base_apex_idx] + DZ[base_apex_idx]
    delta_eff = base_apex_z - sphere_apex_z

    # contact element output (io1): striker = NS7, 625 nodes
    # node IDs in io1 are sequential 1..625 but correspond to global NS7 IDs.
    # exodus stores node_num_map giving global ids
    nmap = f1.variables["node_num_map"][:]   # 1-indexed global ids
    Fz   = f1.variables["vals_nod_var6"][-1]
    Dz1  = f1.variables["vals_nod_var3"][-1]

    # global coordinates of strikers and current radius from axis
    g_cx = cx[nmap - 1]
    g_cy = cy[nmap - 1]
    g_cz = cz[nmap - 1]
    g_DX = DX[nmap - 1]
    g_DY = DY[nmap - 1]
    g_DZ = DZ[nmap - 1]

    cur_x = g_cx + g_DX
    cur_y = g_cy + g_DY
    cur_z = g_cz + g_DZ
    r = np.sqrt(cur_x * cur_x + cur_y * cur_y)

    # F_Z on each striker is the contact force in +z (sphere pushed up by
    # base; n_base_facet = +z, dphi = -fK·h·area > 0 since h<0).  Sum is
    # the upward reaction the base exerts on the sphere quadrant; equals
    # the contact load magnitude.  Full sphere F = 4 × F_quarter.
    F_quarter = np.sum(Fz)

    # contact radius: largest r where Fz exceeds a small fraction of peak
    threshold = 0.05 * Fz.max() if Fz.max() > 0 else 0.0
    in_contact = Fz > threshold
    a_sim = r[in_contact].max() if in_contact.any() else 0.0

    # pressure profile: pressure = +F_z / striker_tributary_area
    striker_area_map, _ = striker_areas_from_geom(geom_path)
    area = np.array([striker_area_map[int(nid)] for nid in nmap])
    p = np.where(area > 0, Fz / area, 0.0)

    return {
        "label": label,
        "delta_eff": delta_eff,
        "F_quarter": F_quarter,
        "a_sim": a_sim,
        "r": r,
        "p": p,
    }


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    geom = os.path.join(here, "geometry", "hertz.geom")

    runs = []
    for stem, label, estar in [("hertz_explicit", "explicit",
                                Estar_for(E, NU, E, NU)),       # both E
                               ("hertz_implicit", "implicit",
                                Estar_for(E, NU, 10 * E, NU))]: # base 10x
        io0 = os.path.join(here, f"{stem}.io0.exo")
        io1 = os.path.join(here, f"{stem}.io1.exo")
        if not os.path.exists(io0) or not os.path.exists(io1):
            print(f"[skip] missing output for {label}: {io0}")
            continue
        d = analyse(io0, io1, geom, label)
        d["Estar"] = estar
        runs.append(d)

    if not runs:
        print("no outputs to compare")
        sys.exit(1)

    # Tabulate vs analytical at the prescribed indenter displacement δ_BC
    print()
    print(f"  R = {R} mm,  ν = {NU},  δ_BC = {DELTA_BC} mm (prescribed)")
    print(f"  {'run':<10} {'E*':>8} {'a_sim':>8} {'a_Hz':>8}"
          f" {'F_quarter_sim':>14} {'F_Hz/4':>10}"
          f" {'p₀_sim':>10} {'p₀_Hz':>10}")
    for r in runs:
        a_h, F_h, p0_h = hertz(DELTA_BC, r["Estar"])
        p0_sim = r["p"].max()
        print(f"  {r['label']:<10} {r['Estar']:>8.0f}"
              f" {r['a_sim']:>8.4f} {a_h:>8.4f}"
              f" {r['F_quarter']:>14.2f} {F_h / 4:>10.2f}"
              f" {p0_sim:>10.1f} {p0_h:>10.1f}")
    print()
    print(f"  (δ_apex_sim = sphere apex DZ - base apex DZ; small because "
          f"both bodies deform.")
    print(f"   Hertz analytical above uses δ_BC = top displacement; the "
          f"contact zone size matches because")
    print(f"   the indentation imposed by a curved indenter pushed down by "
          f"δ produces a = √(R·δ).)")
    print()
    for r in runs:
        print(f"  {r['label']}: δ_apex_sim = {r['delta_eff']*1e3:.2f} μm")

    # plot p(r) overlay vs Hertz at δ_BC
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for r, color, marker in zip(runs, ["C0", "C1"], ["o", "s"]):
            ax.plot(r["r"], r["p"], marker, color=color,
                    label=f"Tahoe {r['label']} (E*={r['Estar']:.0f})",
                    markersize=4, alpha=0.8)
        # Plot analytical for each run's E*
        for r, color in zip(runs, ["C0", "C1"]):
            a_h, F_h, p0_h = hertz(DELTA_BC, r["Estar"])
            r_plot = np.linspace(0, a_h, 200)
            p_plot = p0_h * np.sqrt(np.clip(1.0 - (r_plot / a_h) ** 2, 0, 1))
            ax.plot(r_plot, p_plot, "--", color=color, alpha=0.8, linewidth=2,
                    label=f"Hertz (E*={r['Estar']:.0f}, p₀={p0_h:.0f})")
        ax.set_xlabel("radial distance from contact axis r [mm]")
        ax.set_ylabel("contact pressure p(r) [MPa]")
        ax.set_title(f"Hertz contact: Tahoe vs analytical "
                     f"(R={R} mm, δ={DELTA_BC} mm)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        out = os.path.join(here, "hertz_pressure.png")
        fig.tight_layout()
        fig.savefig(out, dpi=130)
        print(f"\n  pressure profile plot: {out}")
    except ImportError:
        print("matplotlib not installed; skipping plot")


if __name__ == "__main__":
    main()
