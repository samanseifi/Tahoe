#!/usr/bin/env python3
"""Brinell indentation benchmark mesh — #47.

Quarter-symmetry curved-bottom indenter pressed into a flat elastoplastic
base.  Same topology as `level.5/hertz/generate_hertz_mesh.py` (block-
with-spherical-bottom indenter on flat half-space) — see that file for
the design rationale on the box-with-arced-bottom geometry choice
(local apex IS spherical of radius R, which is all Hertz / Brinell
mechanics needs; box mesh avoids degenerate hexes at the apex).

Differences from the Hertz mesh:
  - smaller R and R_box so plasticity dominates over elasticity
  - the indenter ("block 1") will be driven with a near-rigid elastic
    material from XML;  the block ("block 2") gets J2 plasticity with
    cubic-spline hardening (steel-like)
  - 20 graded load steps push the indenter to δ = 0.3 mm ≈ 0.06 R,
    deep into the fully-plastic regime

Geometry (mm):
  R       = 5.0   indenter sphere radius
  R_box   = 2.5   quarter-square edge length
  h_top   = 4.0   indenter top-plane height
  H_base  = 3.0   base block thickness
  gap     = 1e-3  initial contact-pair gap

Node sets:
  1 — indenter top face (z = h_top + gap): driven by δ
  2 — indenter x = 0  face: symmetry, u_x = 0
  3 — indenter y = 0  face: symmetry, u_y = 0
  4 — base bottom face: clamped
  5 — base x = 0  face: symmetry, u_x = 0
  6 — base y = 0  face: symmetry, u_y = 0
  7 — indenter bottom face (contact strikers)
  8 — base top face   (contact strikers, also used for output extraction)

Side sets:
  1 — indenter bottom (block 1, face 1)
  2 — base top        (block 2, face 2)
"""

import math
import os


def grade_tanh(N: int, scale: float = 2.5):
    """N+1 grading values in [0, 1], concentrated at 0 when scale > 0."""
    if scale <= 0:
        return [i / N for i in range(N + 1)]
    pts = [(1.0 - math.tanh(scale * (1.0 - i / N)) / math.tanh(scale))
           for i in range(N + 1)]
    pts[0] = 0.0
    pts[-1] = 1.0
    return pts


def write_mesh(name: str,
               nx: int = 20, ny: int = 20, nz_sph: int = 12,
               mx: int = 20, my: int = 20, mz_base: int = 10,
               R: float = 5.0, R_box: float = 2.5,
               h_top: float = 4.0, H_base: float = 3.0,
               grade_xy: float = 1.5, grade_z: float = 1.5,
               gap: float = 1.0e-3):
    z_bot_corner = R - math.sqrt(max(R * R - 2.0 * R_box * R_box, 0.0))
    if h_top <= z_bot_corner:
        raise ValueError(f"h_top={h_top} must exceed z_bot at corner "
                         f"({z_bot_corner:.4f}); reduce R_box or raise h_top")
    nnx, nny, nnz = nx + 1, ny + 1, nz_sph + 1
    mnx, mny, mnz = mx + 1, my + 1, mz_base + 1

    sph_nnod = nnx * nny * nnz
    base_nnod = mnx * mny * mnz
    numnod = sph_nnod + base_nnod

    sphere_nelem = nx * ny * nz_sph
    base_nelem   = mx * my * mz_base
    numel = sphere_nelem + base_nelem

    outdir = os.path.dirname(os.path.abspath(__file__))
    base = os.path.join(outdir, name)

    xs = [R_box * t for t in grade_tanh(nx, grade_xy)]
    ys = [R_box * t for t in grade_tanh(ny, grade_xy)]
    zs_norm = grade_tanh(nz_sph, grade_z)

    def nid_sphere(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    def nid_base(i, j, k):
        return sph_nnod + k * mnx * mny + j * mnx + i + 1

    coords = []
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                x, y = xs[i], ys[j]
                z_bot = R - math.sqrt(max(R*R - x*x - y*y, 0.0))
                z = z_bot + zs_norm[k] * (h_top - z_bot)
                coords.append((nid_sphere(i, j, k), x, y, z + gap))

    xs_b = [R_box * t for t in grade_tanh(mx, grade_xy)]
    ys_b = [R_box * t for t in grade_tanh(my, grade_xy)]
    for k in range(mnz):
        for j in range(mny):
            for i in range(mnx):
                x, y = xs_b[i], ys_b[j]
                t_top = grade_tanh(mz_base, grade_z)[mz_base - k]
                z = -H_base * t_top
                coords.append((nid_base(i, j, k), x, y, z))

    conn_sphere = []
    eid = 1
    for k in range(nz_sph):
        for j in range(ny):
            for i in range(nx):
                ns = [nid_sphere(i,   j,   k),   nid_sphere(i+1, j,   k),
                      nid_sphere(i+1, j+1, k),   nid_sphere(i,   j+1, k),
                      nid_sphere(i,   j,   k+1), nid_sphere(i+1, j,   k+1),
                      nid_sphere(i+1, j+1, k+1), nid_sphere(i,   j+1, k+1)]
                conn_sphere.append((eid, ns)); eid += 1

    conn_base = []
    eid = 1
    for k in range(mz_base):
        for j in range(my):
            for i in range(mx):
                ns = [nid_base(i,   j,   k),   nid_base(i+1, j,   k),
                      nid_base(i+1, j+1, k),   nid_base(i,   j+1, k),
                      nid_base(i,   j,   k+1), nid_base(i+1, j,   k+1),
                      nid_base(i+1, j+1, k+1), nid_base(i,   j+1, k+1)]
                conn_base.append((eid, ns)); eid += 1

    ns1 = [nid_sphere(i, j, nz_sph) for j in range(nny) for i in range(nnx)]   # indenter top
    ns2 = [nid_sphere(0, j, k)      for k in range(nnz) for j in range(nny)]   # indenter x=0
    ns3 = [nid_sphere(i, 0, k)      for k in range(nnz) for i in range(nnx)]   # indenter y=0
    ns4 = [nid_base(i, j, 0)        for j in range(mny) for i in range(mnx)]   # base bottom
    ns5 = [nid_base(0, j, k)        for k in range(mnz) for j in range(mny)]   # base x=0
    ns6 = [nid_base(i, 0, k)        for k in range(mnz) for i in range(mnx)]   # base y=0
    ns7 = [nid_sphere(i, j, 0)      for j in range(nny) for i in range(nnx)]   # indenter contact
    ns8 = [nid_base(i, j, mz_base)  for j in range(mny) for i in range(mnx)]   # base top contact

    ss1 = []
    eid = 1
    for k in range(nz_sph):
        for j in range(ny):
            for i in range(nx):
                if k == 0: ss1.append((eid, 1))
                eid += 1
    ss2 = []
    eid = 1
    for k in range(mz_base):
        for j in range(my):
            for i in range(mx):
                if k == mz_base - 1: ss2.append((eid, 2))
                eid += 1

    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\nbrinell_indenter_on_plastic_block\n")
        f.write("*dimensions\n")
        f.write(f"{numnod}   # number of nodes\n3   # ndim\n")
        f.write("2   # number of element sets\n")
        f.write(f"1   {sphere_nelem}   8\n")
        f.write(f"2   {base_nelem}   8\n")
        f.write("8   # number of node sets\n")
        for idx, lst in enumerate([ns1, ns2, ns3, ns4, ns5, ns6, ns7, ns8]):
            f.write(f"{idx+1}   {len(lst)}\n")
        f.write("2   # number of side sets\n")
        f.write(f"1  1  {len(ss1)}\n")
        f.write(f"2  2  {len(ss2)}\n")
        f.write("# end dimensions\n")

        f.write("*nodesets\n")
        for lst in [ns1, ns2, ns3, ns4, ns5, ns6, ns7, ns8]:
            f.write("*set\n")
            f.write(f"{len(lst)}\n")
            f.write(" ".join(str(n) for n in lst) + "\n")

        f.write("*sidesets\n")
        for lst in [ss1, ss2]:
            f.write("*set\n")
            f.write(f"{len(lst)}\n")
            for elem, face in lst:
                f.write(f"{elem}  {face}\n")

        f.write("*elements\n")
        for cn in [conn_sphere, conn_base]:
            f.write("*set\n")
            f.write(f"{len(cn)}   # nelem\n")
            f.write("8   # nen\n")
            for eid, ns in cn:
                f.write(f"{eid}  " + "  ".join(str(x) for x in ns) + "\n")

        f.write("*nodes\n")
        f.write(f"{numnod}\n3\n")
        for nid, x, y, z in coords:
            f.write(f"{nid}  {x:.15e}  {y:.15e}  {z:.15e}\n")

    for ext in (".nd", ".es0", ".es1") + tuple(f".ns{i}" for i in range(8)) + (".ss0", ".ss1"):
        p = f"{base}.geom{ext}"
        if os.path.exists(p):
            os.remove(p)

    print(f"{name}: indenter {nx}x{ny}x{nz_sph} = {sphere_nelem} hex,"
          f" block {mx}x{my}x{mz_base} = {base_nelem} hex, total {numel} hex,"
          f" {numnod} nodes")
    print(f"  expected Hertz a_Hz = √(R·δ_yield) at first yield (δ_yield ≈ "
          f"σ_y0·π·R / (2 E*) for σ_y0 = 250 MPa, see README.md)")


if __name__ == "__main__":
    # Smoke variant — fast enough for CI / autonomous loops (~3k hex, ~5 min serial).
    write_mesh("brinell_smoke",
               nx=12, ny=12, nz_sph=8,
               mx=12, my=12, mz_base=8)
    # Fine variant — converges to the Tabor relation; longer runtime (~1 h serial).
    write_mesh("brinell")
