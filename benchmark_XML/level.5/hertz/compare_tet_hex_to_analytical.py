#!/usr/bin/env python3
"""Hertz analytical comparison for the Tet4 indenter / Hex8 base run.

Same physics and analytical formulas as compare_to_analytical.py, but
adapts to the mixed-element-type io layout:
  io0.exo — Tet4 indenter (block 1)
  io1.exo — Hex8 base       (block 2)
  io2.exo — contact strikers (NS7), with F_X / F_Y / F_Z

The original script assumed everything was in one io0 (both blocks
shared an <updated_lagrangian>) plus io1 contact.  Here we just point
the existing analyse() at io0 (sphere field) and io2 (contact strikers).
The geom-parser-driven striker-area computation keeps working unchanged.
"""
import os
import sys

import numpy as np

# Reuse the analytical formulas / fitter / plotter from the original.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_to_analytical as core


# ---------------------------------------------------------------------
# Element-type-aware striker area parser.
#
# core.striker_areas_from_geom only knows Hex8 face numbering.  For the
# Tet4 indenter, SS1 uses face 4 (Tahoe Tet4 1-indexed) which is the
# triangle opposite local node 3, with vertex order {0, 2, 1}.  We
# detect element type per block by `nen` in the .geom *elements section
# and dispatch to the appropriate face_local table.  Tet4 face area is
# the triangle area; nodal share is 1/3 of the triangle area.
# ---------------------------------------------------------------------
def striker_areas_from_geom_mixed(geom_path):
    with open(geom_path) as f:
        text = f.read()
    nodes_idx = text.rfind("\n*nodes\n") + 1
    sidesets_idx = text.index("*sidesets")
    elements_idx = text.index("*elements")

    # parse nodes
    nodes_block = text[nodes_idx:]
    lines = nodes_block.split("\n")
    n_nod = int(lines[1].strip())
    coords = {}
    for i in range(3, 3 + n_nod):
        parts = lines[i].split()
        coords[int(parts[0])] = (float(parts[1]), float(parts[2]),
                                 float(parts[3]))

    # parse element blocks: track nen per block
    elements_block = text[elements_idx:nodes_idx]
    blocks = elements_block.split("*set")[1:]
    elem_conn = {}                # (block_idx, local_eid) -> connectivity
    block_nen = {}                # block_idx -> nen
    for bi, blk in enumerate(blocks):
        ll = blk.strip().split("\n")
        nelem = int(ll[0].split()[0])
        nen = int(ll[1].split()[0])
        block_nen[bi + 1] = nen
        for k in range(nelem):
            parts = ll[2 + k].split()
            eid = int(parts[0])
            elem_conn[(bi + 1, eid)] = list(map(int, parts[1:1 + nen]))

    # parse SS1 (sphere/indenter bottom).  Sideset header gives the
    # element-block ID, but we already know SS1 is on block 1.
    ss_block = text[sidesets_idx:elements_idx]
    ss = ss_block.split("*set")[1:]
    ss1_lines = ss[0].strip().split("\n")
    n_ss1 = int(ss1_lines[0])
    ss1_pairs = []
    for ln in ss1_lines[1:1 + n_ss1]:
        e, fc = ln.split()
        ss1_pairs.append((int(e), int(fc)))

    # face → local-node-list, by element type.
    # Hex8 1-indexed (Tahoe): face 1 = {0,3,2,1}, face 2 = {4,5,6,7}, …
    # Tet4 1-indexed: face 1 = {0,1,3}, face 2 = {1,2,3}, face 3 = {2,0,3}, face 4 = {0,2,1}
    face_local_hex8 = {1: [0, 3, 2, 1], 2: [4, 5, 6, 7]}
    face_local_tet4 = {1: [0, 1, 3], 2: [1, 2, 3],
                       3: [2, 0, 3], 4: [0, 2, 1]}

    striker_area = {}
    nen1 = block_nen[1]
    if nen1 == 8:
        face_local = face_local_hex8
    elif nen1 == 4:
        face_local = face_local_tet4
    else:
        raise ValueError(f"unsupported indenter block nen={nen1}")

    for eid, fc in ss1_pairs:
        ns = elem_conn[(1, eid)]
        face_nodes = [ns[k] for k in face_local[fc]]
        if len(face_nodes) == 4:
            p0 = np.array(coords[face_nodes[0]])
            p1 = np.array(coords[face_nodes[1]])
            p2 = np.array(coords[face_nodes[2]])
            p3 = np.array(coords[face_nodes[3]])
            area = 0.5 * (np.linalg.norm(np.cross(p1 - p0, p2 - p0))
                          + np.linalg.norm(np.cross(p2 - p0, p3 - p0)))
            share = 0.25 * area
        else:                                    # triangle (Tet4 face)
            p0 = np.array(coords[face_nodes[0]])
            p1 = np.array(coords[face_nodes[1]])
            p2 = np.array(coords[face_nodes[2]])
            area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0))
            share = area / 3.0
        for nid in face_nodes:
            striker_area[nid] = striker_area.get(nid, 0.0) + share
    return striker_area, coords


# Monkey-patch the core parser so core.analyse() calls our nen-aware version.
core.striker_areas_from_geom = striker_areas_from_geom_mixed


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    geom = os.path.join(here, "geometry", "hertz_tet_hex.geom")

    runs = []
    for stem, label in [("hertz_tet_hex_explicit", "explicit (Tet+Hex)"),
                        ("hertz_tet_hex_implicit", "implicit (Tet+Hex)")]:
        # Use the merged exo so sphere+base coords are in one place.  The
        # merged file preserves the original global node IDs (sphere=1..n_a,
        # base=n_a+1..n_a+n_b), so io2's node_num_map references still hit
        # the right rows.  Run merge_io.py first if the merged file is
        # missing.
        merged = os.path.join(here, f"{stem}_merged.exo")
        io0 = os.path.join(here, f"{stem}.io0.exo")
        io1 = os.path.join(here, f"{stem}.io1.exo")
        io_contact = os.path.join(here, f"{stem}.io2.exo")
        if not os.path.exists(merged) and os.path.exists(io0) and os.path.exists(io1):
            print(f"  [merge] {os.path.basename(io0)} + "
                  f"{os.path.basename(io1)} → {os.path.basename(merged)}")
            import merge_io
            merge_io.merge(io0, io1, merged)
        if not (os.path.exists(merged) and os.path.exists(io_contact)):
            print(f"[skip] missing output for {label}: {merged} or {io_contact}")
            continue
        d = core.analyse(merged, io_contact, geom, label)
        d["Estar"] = core.ESTAR_TWO
        runs.append(d)

    if not runs:
        print("no outputs to compare")
        sys.exit(1)

    # Hertz analytical at δ = δ_BC (half-space limit reference)
    a_h_BC, F_h_BC, p0_h_BC = core.hertz(core.DELTA_BC, core.ESTAR_TWO)
    p_mid_h_BC = p0_h_BC * (1.0 - 0.5 ** 2) ** 0.5

    print()
    print(f"  R = {core.R} mm,  ν = {core.NU},  δ_BC = {core.DELTA_BC} mm")
    print(f"  E* = E/(2(1-ν²)) = {core.ESTAR_TWO:.0f} MPa  (two identical "
          f"elastic bodies)")
    print()
    print(f"  Hertz @ δ_BC: a = {a_h_BC:.4f} mm,  F/4 = {F_h_BC/4:.1f} N,"
          f"  p₀ = {p0_h_BC:.0f} MPa,  p_mid = {p_mid_h_BC:.1f} MPa")
    print()
    print(f"  Per-run δ_apex (penalty give; small when bodies are half-space-like):")

    for r in runs:
        # Hertz fit
        p0_fit, a_fit, rms = core.fit_hertz(r["r"], r["p"])
        r["p0_fit"] = p0_fit
        r["a_fit"]  = a_fit
        r["rms"]    = rms
        # Mid-radius pressure: 0.4·a_fit to 0.6·a_fit
        import numpy as np
        mask = (r["r"] >= 0.4 * a_fit) & (r["r"] <= 0.6 * a_fit)
        r["p_mid_sim"] = float(r["p"][mask].mean()) if mask.any() else 0.0
        print(f"    {r['label']}: δ_apex = {r['delta_eff']*1e3:.2f} μm")
    print()

    print("  Direct extraction vs Hertz @ δ_BC:")
    print(f"  {'run':<22} {'a_sim':>8} {'F_q_sim':>10} "
          f"{'p_mid_sim':>11} {'p₀_sim':>11}    {'a%':>6} {'F%':>6} "
          f"{'p_mid%':>7} {'p₀%':>6}")
    for r in runs:
        p0_sim = float(r["p"].max())
        ea = (r["a_sim"] - a_h_BC) / a_h_BC * 100
        eF = (r["F_quarter"] - F_h_BC / 4) / (F_h_BC / 4) * 100
        epm = (r["p_mid_sim"] - p_mid_h_BC) / p_mid_h_BC * 100
        ep0 = (p0_sim - p0_h_BC) / p0_h_BC * 100
        print(f"  {r['label']:<22} {r['a_sim']:>8.4f} {r['F_quarter']:>10.1f}"
              f" {r['p_mid_sim']:>11.1f} {p0_sim:>11.1f}    "
              f"{ea:+6.1f} {eF:+6.1f} {epm:+7.1f} {ep0:+6.1f}")
    print()
    print("  Hertz fit (least-squares to p₀√(1−(r/a)²)) vs Hertz @ δ_BC:")
    print(f"  {'run':<22} {'a_fit':>8} {'p₀_fit':>11}  "
          f"{'F_fit/4':>10} {'rms':>10}    {'a%':>6} {'p₀%':>6} {'F%':>6}")
    import math
    for r in runs:
        a_f = r["a_fit"]; p0_f = r["p0_fit"]
        F_f = (2.0 / 3.0) * math.pi * a_f * a_f * p0_f
        ea = (a_f - a_h_BC) / a_h_BC * 100
        ep0 = (p0_f - p0_h_BC) / p0_h_BC * 100
        eF = (F_f / 4 - F_h_BC / 4) / (F_h_BC / 4) * 100
        print(f"  {r['label']:<22} {a_f:>8.4f} {p0_f:>11.1f}  "
              f"{F_f / 4:>10.1f} {r['rms']:>10.1f}    "
              f"{ea:+6.1f} {ep0:+6.1f} {eF:+6.1f}")
    print()

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        ax, ax2 = axes
        for r, color, marker in zip(runs, ["C0", "C1"], ["o", "s"]):
            ax.plot(r["r"], r["p"], marker, color=color,
                    label=f"Tahoe {r['label']}", markersize=4, alpha=0.6)
        r_plot = np.linspace(0, a_h_BC, 200)
        p_plot = p0_h_BC * np.sqrt(np.clip(1.0 - (r_plot / a_h_BC) ** 2, 0, 1))
        ax.plot(r_plot, p_plot, "k--", linewidth=2.5,
                label=f"Hertz analytical (a={a_h_BC:.3f}, p₀={p0_h_BC:.0f})")
        for r, color in zip(runs, ["C0", "C1"]):
            r_plot = np.linspace(0, r["a_fit"], 200)
            p_plot = r["p0_fit"] * np.sqrt(np.clip(
                1.0 - (r_plot / r["a_fit"]) ** 2, 0, 1))
            ax.plot(r_plot, p_plot, "-", color=color, alpha=0.5, linewidth=1.5,
                    label=f"fit {r['label']} (a={r['a_fit']:.3f}, p₀={r['p0_fit']:.0f})")
        ax.set_xlabel("radial distance r [mm]")
        ax.set_ylabel("contact pressure p(r) [MPa]")
        ax.set_title(f"Hertz: Tet4 indenter on Hex8 base "
                     f"(R={core.R} mm, E*={core.ESTAR_TWO:.0f} MPa, δ={core.DELTA_BC} mm)")
        ax.legend()
        ax.grid(True, alpha=0.3)

        for r, color, marker in zip(runs, ["C0", "C1"], ["o", "s"]):
            ax2.plot(r["r"] / a_h_BC, r["p"] / p0_h_BC, marker, color=color,
                     label=f"Tahoe {r['label']}", markersize=4, alpha=0.7)
        x = np.linspace(0, 1, 200)
        ax2.plot(x, np.sqrt(np.clip(1 - x ** 2, 0, 1)), "k--", linewidth=2,
                 label="Hertz analytical")
        ax2.set_xlabel("r / a_Hz"); ax2.set_ylabel("p / p₀_Hz")
        ax2.set_xlim(0, 2.0); ax2.set_ylim(0, 2.0)
        ax2.axhline(1.0, color="grey", lw=0.6, ls=":")
        ax2.axvline(1.0, color="grey", lw=0.6, ls=":")
        ax2.set_title("Normalised: collapses onto universal Hertz curve")
        ax2.legend(); ax2.grid(True, alpha=0.3)

        out = os.path.join(here, "hertz_tet_hex_pressure.png")
        fig.tight_layout()
        fig.savefig(out, dpi=130)
        print(f"  pressure profile plot: {out}")
    except ImportError:
        print("matplotlib not installed; skipping plot")


if __name__ == "__main__":
    main()
