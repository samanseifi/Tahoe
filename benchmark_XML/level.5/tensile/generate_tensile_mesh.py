#!/usr/bin/env python3
"""Tensile-test mesh — quarter-symmetry bar (#49).

Pulls a slim rectangular bar in tension to verify that <Simo_J2> with a
cubic-spline isotropic hardening curve reproduces the input σ_y(ε_p)
data under uniaxial loading.

Geometry (mm), quarter-symmetry:
  L_x = 0.5    (half-width)
  L_y = 0.5    (half-depth)
  L_z = 5.0    (half-length; full bar would be 10 mm)

Mesh:
  nx = ny = 4,  nz = 40   ⇒ 640 Hex8

Node sets:
  1 — top face (z = L_z)     prescribed u_z = δ
  2 — bottom face (z = 0)     u_z = 0
  3 — x = 0 face              symmetry, u_x = 0
  4 — y = 0 face              symmetry, u_y = 0
  5 — single corner node at (L_x, L_y, L_z/2)   used by post-processor
       to pick a deterministic IP for σ-α extraction

The output (io0.exo) carries D_X/Y/Z, Cauchy stress s11..s12, plus the
J2 material outputs alpha, norm_beta, VM_Kirch, press.
"""

import os


def write_mesh(name: str = "tensile",
               nx: int = 4, ny: int = 4, nz: int = 40,
               Lx: float = 0.5, Ly: float = 0.5, Lz: float = 5.0):
    outdir = os.path.dirname(os.path.abspath(__file__))
    base = os.path.join(outdir, name)

    nnx, nny, nnz = nx + 1, ny + 1, nz + 1
    numnod = nnx * nny * nnz
    numel  = nx * ny * nz

    def nid(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    # coordinates: uniform grid
    coords = []
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                x = Lx * i / nx
                y = Ly * j / ny
                z = Lz * k / nz
                coords.append((nid(i, j, k), x, y, z))

    # element connectivity
    conn = []
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ns = [nid(i,   j,   k),   nid(i+1, j,   k),
                      nid(i+1, j+1, k),   nid(i,   j+1, k),
                      nid(i,   j,   k+1), nid(i+1, j,   k+1),
                      nid(i+1, j+1, k+1), nid(i,   j+1, k+1)]
                conn.append((eid, ns)); eid += 1

    # node sets
    ns1 = [nid(i, j, nz) for j in range(nny) for i in range(nnx)]   # top
    ns2 = [nid(i, j, 0)  for j in range(nny) for i in range(nnx)]   # bottom
    ns3 = [nid(0, j, k)  for k in range(nnz) for j in range(nny)]   # x=0
    ns4 = [nid(i, 0, k)  for k in range(nnz) for i in range(nnx)]   # y=0
    # NS5 — a deterministic corner-strip used as the σ–α probe by the post-
    # processor.  We pick the centre-column (i=nx, j=ny) so the strip stays
    # uniaxial throughout the deformation (centroidal axis is symmetric).
    ns5 = [nid(nx, ny, k) for k in range(nnz)]

    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\ntensile_quarter_bar\n")
        f.write("*dimensions\n")
        f.write(f"{numnod}   # number of nodes\n3   # ndim\n")
        f.write("1   # number of element sets\n")
        f.write(f"1   {numel}   8\n")
        f.write("5   # number of node sets\n")
        for idx, lst in enumerate([ns1, ns2, ns3, ns4, ns5]):
            f.write(f"{idx+1}   {len(lst)}\n")
        f.write("0   # number of side sets\n")
        f.write("# end dimensions\n")

        f.write("*nodesets\n")
        for lst in [ns1, ns2, ns3, ns4, ns5]:
            f.write("*set\n")
            f.write(f"{len(lst)}\n")
            f.write(" ".join(str(n) for n in lst) + "\n")

        f.write("*sidesets\n")

        f.write("*elements\n*set\n")
        f.write(f"{numel}   # nelem\n8   # nen\n")
        for eid, ns in conn:
            f.write(f"{eid}  " + "  ".join(str(x) for x in ns) + "\n")

        f.write("*nodes\n")
        f.write(f"{numnod}\n3\n")
        for n, x, y, z in coords:
            f.write(f"{n}  {x:.15e}  {y:.15e}  {z:.15e}\n")

    # clean stale per-set files
    for ext in (".nd", ".es0") + tuple(f".ns{i}" for i in range(5)):
        p = f"{base}.geom{ext}"
        if os.path.exists(p):
            os.remove(p)

    print(f"{name}: {nx}x{ny}x{nz} = {numel} Hex8, {numnod} nodes")


if __name__ == "__main__":
    write_mesh()
