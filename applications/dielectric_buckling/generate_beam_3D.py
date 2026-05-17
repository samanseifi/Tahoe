#!/usr/bin/env python3
"""
Generate a perfect 3D hex8 beam mesh (no geometric imperfection).

Domain: x in [-10,10], y in [-0.5,0.5], z in [0,Lz]
  Nx=40 elements in x  (dx = 0.5)
  Ny=2  elements in y  (dy = 0.5, thickness / electric-field direction)
  Nz    elements in z  (width, default 2 → square cross-section)

Electric field applied across y: bottom face (y=-0.5) = ground,
                                  top face   (y=+0.5) = driven electrode.
Buckling is seeded by floating-point round-off in the time integrator,
not by a geometric imperfection, so the mesh is perfectly regular.

Node numbering: k*(Nx+1)*(Ny+1) + j*(Nx+1) + i + 1
  i: x-index (0..Nx), j: y-index (0..Ny), k: z-index (0..Nz)

Node sets:
  1  y=-0.5   ground electrode      (Psi = 0)
  2  x=+10    right clamp           (D_X=D_Y=D_Z=0)
  3  y=+0.5   driven electrode      (Psi = V)
  4  x=-10    left clamp            (D_X=D_Y=D_Z=0)
  5  z=0      side roller           (D_Z = 0)
  6  z=Lz     side roller           (D_Z = 0)
"""

import argparse
import os


def generate(Lx, Lz, Nx, Ny, Nz, out_path):
    Ly = 1.0          # thickness in y — fixed to match 2D
    x0 = -Lx / 2.0

    nnx, nny, nnz = Nx + 1, Ny + 1, Nz + 1
    num_nodes = nnx * nny * nnz
    num_elem  = Nx  * Ny  * Nz

    def nid(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    # Nodes — perfect geometry; numerical round-off in the HHT integrator
    # seeds the bifurcation without baking in a preferred mode shape.
    nodes = []
    for k in range(nnz):
        z = k * Lz / Nz
        for j in range(nny):
            y = -Ly / 2.0 + j * Ly / Ny
            for i in range(nnx):
                x = x0 + i * Lx / Nx
                nodes.append((nid(i, j, k), x, y, z))

    # Elements  (connectivity matches Tahoe hex8 single_hex_3D.geom pattern)
    elements = []
    eid = 1
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                conn = (
                    nid(i,   j,   k  ), nid(i+1, j,   k  ),
                    nid(i+1, j+1, k  ), nid(i,   j+1, k  ),
                    nid(i,   j,   k+1), nid(i+1, j,   k+1),
                    nid(i+1, j+1, k+1), nid(i,   j+1, k+1),
                )
                elements.append((eid, *conn))
                eid += 1

    # Node sets
    ns = {sid: [] for sid in range(1, 7)}
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                n = nid(i, j, k)
                if j == 0:    ns[1].append(n)
                if i == Nx:   ns[2].append(n)
                if j == Ny:   ns[3].append(n)
                if i == 0:    ns[4].append(n)
                if k == 0:    ns[5].append(n)
                if k == Nz:   ns[6].append(n)

    ns_labels = {
        1: "y=-0.5  ground electrode (Psi=0)",
        2: "x=+10   right clamp",
        3: "y=+0.5  driven electrode (Psi=V)",
        4: "x=-10   left clamp",
        5: "z=0     side roller",
        6: f"z={Lz}  side roller",
    }

    with open(out_path, "w") as f:
        f.write("*version\n1.0\n")
        f.write("*title\n")
        f.write(f"3D beam: {Nx}x{Ny}x{Nz} hex8, Lx={Lx} Ly={Ly} Lz={Lz}\n")
        f.write("*dimensions\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("3   # number of spatial dimensions\n")
        f.write("1   # number of element sets\n")
        f.write("# [ID] [nel] [nen]\n")
        f.write(f"1  {num_elem}  8\n")
        f.write(f"{len(ns)}   # number of node sets\n")
        f.write("# [ID] [nnd]\n")
        for sid in sorted(ns):
            f.write(f"{sid}  {len(ns[sid])}\n")
        f.write("0  # number of side sets\n")
        f.write("# end dimensions\n")

        f.write("*nodesets\n")
        for sid in sorted(ns):
            lst = ns[sid]
            f.write("*set\n")
            f.write(f"{len(lst)}  # number of nodes\n")
            for k0 in range(0, len(lst), 10):
                f.write("  ".join(str(n) for n in lst[k0:k0 + 10]) + "\n")
        f.write("# end node sets\n")

        f.write("*sidesets\n# end side sets\n")

        f.write("*elements\n*set\n")
        f.write(f"{num_elem}  # number of elements\n")
        f.write("8  # number of element nodes\n")
        for e in elements:
            f.write("  " + "  ".join(str(v) for v in e) + "\n")
        f.write("# end elements\n")

        f.write("*nodes\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("3   # number of spatial dimensions\n")
        for nd in nodes:
            f.write(f"  {nd[0]}  {nd[1]:.10e}  {nd[2]:.10e}  {nd[3]:.10e}\n")

    print(f"Written: {out_path}")
    print(f"  Nodes: {num_nodes}   Elements: {num_elem}")
    for sid in sorted(ns):
        print(f"  NS {sid} ({ns_labels[sid]}): {len(ns[sid])} nodes")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lx",  type=float, default=20.0,
                    help="beam length (default 20)")
    ap.add_argument("--Lz",  type=float, default=1.0,
                    help="beam width in z (default 1.0 → square cross-section)")
    ap.add_argument("--Nx",  type=int,   default=40,
                    help="elements along x (default 40)")
    ap.add_argument("--Ny",  type=int,   default=2,
                    help="elements through y-thickness (default 2)")
    ap.add_argument("--Nz",  type=int,   default=2,
                    help="elements along z-width (default 2)")
    ap.add_argument("--out", default="beam_3D.geom")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)
    generate(args.Lx, args.Lz, args.Nx, args.Ny, args.Nz, out)
