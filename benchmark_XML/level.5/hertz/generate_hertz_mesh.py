#!/usr/bin/env python3
"""Hertz contact benchmark mesh.

Quarter-symmetry hemisphere (deformable) sitting on a stiff base block.
Both meshed as Hex8.  Inline Tahoe geom format with two element blocks
(block 1 = sphere, block 2 = base) and side sets for contact.

Geometry (mm):
  R       = 10.0   sphere radius
  R_box   = 0.5*R  quarter-square edge length (stays well clear of the
                   sphere boundary so element thicknesses stay healthy;
                   the far corner has z_bot = R - sqrt(R² - R_box²·2)
                   ≈ 0.29R, giving slab thickness ~0.71R there)
  H_top   = R      sphere apex-to-top distance (top is flat at z=R)

Sphere geometry (block 1):
  - (x, y) ∈ [0, R_box] × [0, R_box]    quarter-square
  - z(x, y) bottom = R - sqrt(R² - x² - y²)   (curved, on the sphere surface)
  - z       top    = R                          (flat)

Base block (block 2):
  - (x, y) ∈ [0, R_box] × [0, R_box]
  - z ∈ [-H_base, 0]    flat top at z=0, fully clamped at z=-H_base

Mesh grading (tanh-based) concentrates resolution at:
  - sphere bottom (the contact apex at x=y=0)
  - base top (the contact face)

Node-set IDs:
  1 — sphere top face (z=R): driven downward by prescribed δ
  2 — sphere x=0 face: symmetry (u_x=0)
  3 — sphere y=0 face: symmetry (u_y=0)
  4 — base bottom face (z=-H_base): clamped
  5 — base x=0 face: symmetry (u_x=0)
  6 — base y=0 face: symmetry (u_y=0)
  7 — sphere bottom (contact-side, on sphere surface)
  8 — base top (contact-side, z=0)

Side-set IDs (for contact):
  1 — sphere bottom face   (block 1, face 1 = -z in Tahoe Hex8 1-indexed)
  2 — base top face        (block 2, face 2 = +z in Tahoe Hex8 1-indexed)
"""
import os
import math


def grade_tanh(N: int, scale: float = 2.5):
    """Return N+1 grading values in [0, 1] concentrated at 0.

    Uses h(t) = 1 - tanh(scale*(1-t))/tanh(scale): fine spacing at t=0,
    coarse at t=1.  scale > 0: stronger refinement.  scale=0: uniform.
    """
    if scale <= 0:
        return [i / N for i in range(N + 1)]
    pts = [(1.0 - math.tanh(scale * (1.0 - i / N)) / math.tanh(scale))
           for i in range(N + 1)]
    pts[0] = 0.0
    pts[-1] = 1.0
    return pts


def write_mesh(name: str,
               nx: int = 24, ny: int = 24, nz_sph: int = 12,
               mx: int = 24, my: int = 24, mz_base: int = 6,
               R: float = 10.0, H_base: float = 5.0,
               grade_xy: float = 2.5, grade_z: float = 2.5,
               gap: float = 1.0e-3):
    """gap > 0 lifts the sphere block by `gap` mm so the initial
       configuration has no penetration anywhere on the contact pair —
       penalty contact activates only as compression begins."""
    R_box = 0.5 * R
    nnx, nny, nnz = nx + 1, ny + 1, nz_sph + 1
    mnx, mny, mnz = mx + 1, my + 1, mz_base + 1

    cube1_nnod = nnx * nny * nnz
    cube2_nnod = mnx * mny * mnz
    numnod = cube1_nnod + cube2_nnod

    sphere_nelem = nx * ny * nz_sph
    base_nelem   = mx * my * mz_base
    numel = sphere_nelem + base_nelem

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)

    # ---- sphere coords (block 1) ----
    # graded x,y in [0, R_box]; graded z in [z_bot(x,y), R]
    xs = [R_box * t for t in grade_tanh(nx, grade_xy)]
    ys = [R_box * t for t in grade_tanh(ny, grade_xy)]
    zs_norm = grade_tanh(nz_sph, grade_z)   # normalised, 0 = bottom

    def nid_sphere(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    def nid_base(i, j, k):
        return cube1_nnod + k * mnx * mny + j * mnx + i + 1

    coords = []
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                x, y = xs[i], ys[j]
                z_bot = R - math.sqrt(max(R*R - x*x - y*y, 0.0))
                z_top = R
                z = z_bot + zs_norm[k] * (z_top - z_bot)
                # lift entire sphere by `gap` so initial contact pair
                # has a uniform tiny gap rather than zero penetration
                coords.append((nid_sphere(i, j, k), x, y, z + gap))

    # ---- base coords (block 2) ----
    # uniform x,y; z ∈ [-H_base, 0] graded fine near top
    xs_b = [R_box * t for t in grade_tanh(mx, grade_xy)]
    ys_b = [R_box * t for t in grade_tanh(my, grade_xy)]
    zs_b_norm = grade_tanh(mz_base, grade_z)   # 0 = bottom (clamped)
    for k in range(mnz):
        for j in range(mny):
            for i in range(mnx):
                x, y = xs_b[i], ys_b[j]
                # we want z from -H_base (k=0) to 0 (k=mz_base)
                # finer near z=0 → use 1-tanh from top
                t_top = grade_tanh(mz_base, grade_z)[mz_base - k]   # reverse
                z = -H_base * t_top
                coords.append((nid_base(i, j, k), x, y, z))

    # ---- element connectivity ----
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

    # ---- node sets ----
    ns1 = [nid_sphere(i, j, nz_sph) for j in range(nny) for i in range(nnx)]   # sphere top
    ns2 = [nid_sphere(0, j, k)     for k in range(nnz) for j in range(nny)]   # sphere x=0
    ns3 = [nid_sphere(i, 0, k)     for k in range(nnz) for i in range(nnx)]   # sphere y=0
    ns4 = [nid_base(i, j, 0)        for j in range(mny) for i in range(mnx)]   # base bottom
    ns5 = [nid_base(0, j, k)        for k in range(mnz) for j in range(mny)]   # base x=0
    ns6 = [nid_base(i, 0, k)        for k in range(mnz) for i in range(mnx)]   # base y=0
    ns7 = [nid_sphere(i, j, 0)     for j in range(nny) for i in range(nnx)]   # sphere contact
    ns8 = [nid_base(i, j, mz_base)  for j in range(mny) for i in range(mnx)]   # base contact

    # ---- side sets (block-local element IDs, 1-indexed face IDs) ----
    ss1 = []  # sphere bottom — face 1 (-z), elements at k=0
    eid = 1
    for k in range(nz_sph):
        for j in range(ny):
            for i in range(nx):
                if k == 0: ss1.append((eid, 1))
                eid += 1
    ss2 = []  # base top — face 2 (+z), elements at k=mz_base-1
    eid = 1
    for k in range(mz_base):
        for j in range(my):
            for i in range(mx):
                if k == mz_base - 1: ss2.append((eid, 2))
                eid += 1

    # ---- write inline geom ----
    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\nhertz_quarter_sphere_on_stiff_base\n")
        f.write("*dimensions\n")
        f.write(f"{numnod}   # number of nodes\n3   # ndim\n")
        f.write("2   # number of element sets\n")
        f.write(f"1   {sphere_nelem}   8\n")
        f.write(f"2   {base_nelem}   8\n")
        f.write(f"8   # number of node sets\n")
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

    # cleanup external files (in case of stale)
    for ext in (".nd", ".es0", ".es1") + tuple(f".ns{i}" for i in range(8)) + (".ss0", ".ss1"):
        p = f"{base}.geom{ext}"
        if os.path.exists(p):
            os.remove(p)

    print(f"{name}: sphere {nx}x{ny}x{nz_sph} = {sphere_nelem} hex,"
          f" base {mx}x{my}x{mz_base} = {base_nelem} hex, total {numel} hex,"
          f" {numnod} nodes")


if __name__ == "__main__":
    # Mild grading (scale=1.0): smallest dx ~0.14 mm at apex, dx ~0.55 mm
    # at far corner; ~10 elements across the analytical contact radius
    # a≈0.7 mm.  Keeps Δt ~0.025 μs (5e-3 of total ramp T=100 μs).
    write_mesh("hertz",
               nx=24, ny=24, nz_sph=12,
               mx=24, my=24, mz_base=8,
               grade_xy=1.0, grade_z=1.0)
