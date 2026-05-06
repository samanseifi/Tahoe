#!/usr/bin/env python3
"""Hertz contact benchmark — Tet4 indenter on Hex8 base.

Mirrors the geometry of generate_hertz_mesh.py (curved-bottom block,
quarter symmetry, h_top = H_base = 6 mm so each body is in the half-
space limit Hertz assumes), but the indenter (block 1) is meshed as
Tet4 via the standard 5-tet split of each hex cell.  The base (block 2)
stays Hex8.

Why this is interesting:
  - Demonstrates mixed-element contact: Tet4 facets (3-node triangles)
    pairing with Hex8 facets (4-node quads → split to 2 triangles by
    Contact3DT::ConvertQuadToTri).  The contact element doesn't care
    about bulk element type; it works on side-sets.
  - Drives the BonetTet (ANP-Tet4 / LS-DYNA ELFORM=13) implicit element
    on a real contact problem.  Plain Tet4 locks badly under volumetric
    constraint; F-bar should give Hertz-pressure profile match.

Geometry (mm) — same as the all-hex Hertz benchmark:
  R       = 10.0   sphere radius (defines the curved bottom)
  R_box   = 3.0    quarter-square cross-section edge
  h_top   = 6.0    indenter top-plane height (≈ 8.5 × a_Hz)
  H_base  = 6.0    base block thickness     (≈ 8.5 × a_Hz)

5-tet split of each hex cell (Tahoe-convention corner numbering 0..7,
bottom 0,1,2,3 then top 4,5,6,7):
    tet 0 = (0, 1, 2, 5)   ← contains hex bottom-face nodes 0,1,2
    tet 1 = (0, 2, 3, 7)   ← contains hex bottom-face nodes 0,2,3
    tet 2 = (0, 5, 7, 4)
    tet 3 = (2, 5, 6, 7)
    tet 4 = (0, 2, 7, 5)   central
Tets 0 and 1 own the two triangles that tile the original hex bottom
quad.  Tahoe Tet4 face-4 (1-indexed) = local nodes {0,2,1} — the face
opposite local-3 — has outward normal in the −z direction for a
positively-oriented tet of this layout, which is exactly the direction
we need for the contact normal at the indenter bottom.

Node-set IDs (same as the all-hex benchmark, identical convention):
  1 — sphere top face (z=h_top): driven by prescribed δ
  2 — sphere x=0 face: symmetry (u_x=0)
  3 — sphere y=0 face: symmetry (u_y=0)
  4 — base bottom face (z=−H_base): clamped
  5 — base x=0 face: symmetry
  6 — base y=0 face: symmetry
  7 — sphere bottom (contact strikers, on sphere surface)
  8 — base top (contact face nodes, free)

Side-set IDs:
  1 — sphere bottom triangles  (block 1, 1-indexed face 4 — Tet4)
  2 — base top quads           (block 2, 1-indexed face 2 — Hex8)
"""
import os
import math


# Standard 5-tet split of a hex8 (matches generate_tet_mesh.py).  Hex local
# corner ordering: bottom 0..3, top 4..7.
HEX_TO_TETS = [
    (0, 1, 2, 5),
    (0, 2, 3, 7),
    (0, 5, 7, 4),
    (2, 5, 6, 7),
    (0, 2, 7, 5),  # central — note 7,5 ordering for positive Jacobian
]

# For each tet in the split, list the local-tet face index (0-based) whose
# triangle covers part of the hex's bottom face (the k=0 face of the cell).
# −1 means "this tet does not touch the bottom face."
# Tahoe face numbering (TetrahedronT::NodesOnFacet): facet 3 = nodes {0,2,1}.
# Verified: outward normal for face 3 of a positively-oriented tet is the
# average of the cross-products formed by the bottom hex layout, pointing
# in −z for tets 0 and 1.
TET_BOTTOM_FACE = [3, 3, -1, -1, -1]   # 0-indexed; +1 written into geom


def grade_tanh(N: int, scale: float = 1.5):
    """Grading values in [0, 1] concentrated at 0.  See generate_hertz_mesh.py."""
    if scale <= 0:
        return [i / N for i in range(N + 1)]
    pts = [(1.0 - math.tanh(scale * (1.0 - i / N)) / math.tanh(scale))
           for i in range(N + 1)]
    pts[0] = 0.0
    pts[-1] = 1.0
    return pts


def write_mesh(name: str = "hertz_tet_hex",
               nx: int = 24, ny: int = 24, nz_sph: int = 16,
               mx: int = 24, my: int = 24, mz_base: int = 12,
               R: float = 10.0, R_box: float = 3.0,
               h_top: float = 6.0, H_base: float = 6.0,
               grade_xy: float = 1.5, grade_z: float = 1.5,
               gap: float = 1.0e-3):
    """Same xy/z grid as generate_hertz_mesh.py — node coordinates are
    identical so a straight comparison against the all-hex Hertz run is
    valid."""
    z_bot_corner = R - math.sqrt(max(R * R - 2.0 * R_box * R_box, 0.0))
    if h_top <= z_bot_corner:
        raise ValueError(f"h_top={h_top} must exceed z_bot at corner "
                         f"({z_bot_corner:.4f})")

    nnx, nny, nnz_sph = nx + 1, ny + 1, nz_sph + 1
    mnx, mny, mnz_b   = mx + 1, my + 1, mz_base + 1
    cube1_nnod = nnx * nny * nnz_sph
    cube2_nnod = mnx * mny * mnz_b
    numnod = cube1_nnod + cube2_nnod

    sphere_nelem_hex = nx * ny * nz_sph
    sphere_nelem_tet = sphere_nelem_hex * len(HEX_TO_TETS)
    base_nelem_hex   = mx * my * mz_base
    numel = sphere_nelem_tet + base_nelem_hex

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry")
    os.makedirs(outdir, exist_ok=True)
    base_path = os.path.join(outdir, name)

    # ---- coordinates: sphere block (curved bottom + flat top, with gap) ----
    xs = [R_box * t for t in grade_tanh(nx, grade_xy)]
    ys = [R_box * t for t in grade_tanh(ny, grade_xy)]
    zs_norm = grade_tanh(nz_sph, grade_z)

    def nid_sphere(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    def nid_base(i, j, k):
        return cube1_nnod + k * mnx * mny + j * mnx + i + 1

    coords = []
    for k in range(nnz_sph):
        for j in range(nny):
            for i in range(nnx):
                x, y = xs[i], ys[j]
                z_bot = R - math.sqrt(max(R*R - x*x - y*y, 0.0))
                z = z_bot + zs_norm[k] * (h_top - z_bot)
                coords.append((nid_sphere(i, j, k), x, y, z + gap))

    # ---- coordinates: base block (flat slab, clamped at z=−H_base) ----
    xs_b = [R_box * t for t in grade_tanh(mx, grade_xy)]
    ys_b = [R_box * t for t in grade_tanh(my, grade_xy)]
    for k in range(mnz_b):
        for j in range(mny):
            for i in range(mnx):
                t_top = grade_tanh(mz_base, grade_z)[mz_base - k]
                z = -H_base * t_top
                coords.append((nid_base(i, j, k), xs_b[i], ys_b[j], z))

    # ---- connectivity: tet split of sphere block ----
    # Each hex cell at (i,j,k) becomes 5 tets; their global tet IDs follow
    # the order [hex-major, then tet-minor] so tet 0..4 of cell (i,j,k)
    # have IDs hex_id * 5 + 0..4 with hex_id = k*nx*ny + j*nx + i.
    conn_sphere_tet = []
    eid = 1
    for k in range(nz_sph):
        for j in range(ny):
            for i in range(nx):
                hex_corners = [
                    nid_sphere(i,   j,   k),
                    nid_sphere(i+1, j,   k),
                    nid_sphere(i+1, j+1, k),
                    nid_sphere(i,   j+1, k),
                    nid_sphere(i,   j,   k+1),
                    nid_sphere(i+1, j,   k+1),
                    nid_sphere(i+1, j+1, k+1),
                    nid_sphere(i,   j+1, k+1),
                ]
                for tet_local in HEX_TO_TETS:
                    ns = [hex_corners[c] for c in tet_local]
                    conn_sphere_tet.append((eid, ns))
                    eid += 1

    # ---- connectivity: hex base ----
    conn_base = []
    eid = 1
    for k in range(mz_base):
        for j in range(my):
            for i in range(mx):
                ns = [
                    nid_base(i,   j,   k),   nid_base(i+1, j,   k),
                    nid_base(i+1, j+1, k),   nid_base(i,   j+1, k),
                    nid_base(i,   j,   k+1), nid_base(i+1, j,   k+1),
                    nid_base(i+1, j+1, k+1), nid_base(i,   j+1, k+1),
                ]
                conn_base.append((eid, ns))
                eid += 1

    # ---- node sets (same convention as the all-hex Hertz mesh) ----
    ns1 = [nid_sphere(i, j, nz_sph) for j in range(nny) for i in range(nnx)]
    ns2 = [nid_sphere(0, j, k)      for k in range(nnz_sph) for j in range(nny)]
    ns3 = [nid_sphere(i, 0, k)      for k in range(nnz_sph) for i in range(nnx)]
    ns4 = [nid_base(i, j, 0)         for j in range(mny) for i in range(mnx)]
    ns5 = [nid_base(0, j, k)         for k in range(mnz_b) for j in range(mny)]
    ns6 = [nid_base(i, 0, k)         for k in range(mnz_b) for i in range(mnx)]
    ns7 = [nid_sphere(i, j, 0)      for j in range(nny) for i in range(nnx)]
    ns8 = [nid_base(i, j, mz_base)   for j in range(mny) for i in range(mnx)]

    # ---- side sets ----
    # SS1 = sphere bottom triangles.  For each (i,j) at k=0, tets 0 and 1
    # of the cell own the two triangles that tile that hex bottom.  Block-
    # local element IDs in the sphere block: hex_id * 5 + tet_local.
    ss1 = []
    for j in range(ny):
        for i in range(nx):
            hex_id_local = (j * nx + i)              # k=0 layer
            for split_idx, face_local in enumerate(TET_BOTTOM_FACE):
                if face_local < 0:
                    continue
                local_eid = hex_id_local * len(HEX_TO_TETS) + split_idx + 1
                ss1.append((local_eid, face_local + 1))   # 1-indexed face

    # SS2 = base top quads.  Block-local element IDs in the base block.
    ss2 = []
    eid = 1
    for k in range(mz_base):
        for j in range(my):
            for i in range(mx):
                if k == mz_base - 1:
                    ss2.append((eid, 2))
                eid += 1

    # ---- write .geom ----
    with open(f"{base_path}.geom", "w") as f:
        f.write("*version\n1.0\n*title\nhertz_tet_indenter_on_hex_base\n")
        f.write("*dimensions\n")
        f.write(f"{numnod}   # number of nodes\n3   # ndim\n")
        f.write("2   # number of element sets\n")
        f.write(f"1   {sphere_nelem_tet}   4    # block 1: Tet4 indenter\n")
        f.write(f"2   {base_nelem_hex}   8    # block 2: Hex8 base\n")
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
        # block 1: Tet4
        f.write("*set\n")
        f.write(f"{len(conn_sphere_tet)}   # nelem\n")
        f.write("4   # nen (Tet4)\n")
        for eid, ns in conn_sphere_tet:
            f.write(f"{eid}  " + "  ".join(str(x) for x in ns) + "\n")
        # block 2: Hex8
        f.write("*set\n")
        f.write(f"{len(conn_base)}   # nelem\n")
        f.write("8   # nen (Hex8)\n")
        for eid, ns in conn_base:
            f.write(f"{eid}  " + "  ".join(str(x) for x in ns) + "\n")

        f.write("*nodes\n")
        f.write(f"{numnod}\n3\n")
        for nid, x, y, z in coords:
            f.write(f"{nid}  {x:.15e}  {y:.15e}  {z:.15e}\n")

    # cleanup any stale external sub-files (we use inline format)
    for ext in (".nd", ".es0", ".es1") + tuple(f".ns{i}" for i in range(8)) + (".ss0", ".ss1"):
        p = f"{base_path}.geom{ext}"
        if os.path.exists(p):
            os.remove(p)

    print(f"{name}: sphere {nx}x{ny}x{nz_sph}=({sphere_nelem_hex} hex) split into "
          f"{sphere_nelem_tet} tet, base {mx}x{my}x{mz_base}={base_nelem_hex} hex, "
          f"total {numel} elements, {numnod} nodes")
    print(f"   SS1: {len(ss1)} tet triangles, SS2: {len(ss2)} hex quads")


if __name__ == "__main__":
    write_mesh("hertz_tet_hex")
