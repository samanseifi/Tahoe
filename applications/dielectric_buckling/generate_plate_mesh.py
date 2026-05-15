#!/usr/bin/env python3
"""Generate a 3D hex8 plate mesh for dielectric-elastomer buckling.

Domain (default): x in [0,Lx], y in [0,Ly], z in [0,Lz]
  Lx = 20  long axis (slender direction)
  Ly =  2  width
  Lz =  1  thickness  -- electric field is applied across this dimension

A small mode-1 sinusoidal Z imperfection is baked into the coordinates
to seed buckling: z_new = z + amp * sin(pi * x / Lx). It vanishes at the
clamped ends (x=0, x=Lx).

Node sets:
  1  x = 0      clamped end (fix D_X, D_Y, D_Z; tie Psi to ground)
  2  x = Lx     clamped end (fix D_X, D_Y, D_Z)
  3  y = 0      side roller  (fix D_Y)
  4  y = Ly     side roller  (fix D_Y)
  5  z = 0      bottom electrode (Psi = 0)
  6  z = Lz     top electrode    (Psi = V(t))

Output:
  plate.geom  -- Tahoe single-file ASCII geometry
"""
import argparse
import math
import os


def generate(Lx, Ly, Lz, nx, ny, nz, imperfection, out_path):
    nnx, nny, nnz = nx + 1, ny + 1, nz + 1
    num_nodes = nnx * nny * nnz
    num_elem = nx * ny * nz

    def nid(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    nodes = []
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                x = i * Lx / nx
                y = j * Ly / ny
                z = k * Lz / nz
                z += imperfection * math.sin(math.pi * x / Lx)
                nodes.append((nid(i, j, k), x, y, z))

    elements = []
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n1 = nid(i,     j,     k)
                n2 = nid(i + 1, j,     k)
                n3 = nid(i + 1, j + 1, k)
                n4 = nid(i,     j + 1, k)
                n5 = nid(i,     j,     k + 1)
                n6 = nid(i + 1, j,     k + 1)
                n7 = nid(i + 1, j + 1, k + 1)
                n8 = nid(i,     j + 1, k + 1)
                elements.append((eid, n1, n2, n3, n4, n5, n6, n7, n8))
                eid += 1

    ns = {1: [], 2: [], 3: [], 4: [], 5: [], 6: []}
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                n = nid(i, j, k)
                if i == 0:  ns[1].append(n)
                if i == nx: ns[2].append(n)
                if j == 0:  ns[3].append(n)
                if j == ny: ns[4].append(n)
                if k == 0:  ns[5].append(n)
                if k == nz: ns[6].append(n)

    with open(out_path, "w") as f:
        f.write("*version\n1.0\n")
        f.write("*title\n")
        f.write(
            f"Dielectric elastomer plate: {nx}x{ny}x{nz} hex8, "
            f"Lx={Lx} Ly={Ly} Lz={Lz}, imperfection={imperfection}\n"
        )

        f.write("*dimensions\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("3   # number of spatial dimensions\n")
        f.write("1   # number of element sets\n")
        f.write("# [ID] [nel] [nen]\n")
        f.write(f"1   {num_elem}   8\n")
        f.write(f"{len(ns)}   # number of node sets\n")
        f.write("# [ID] [nnd]\n")
        for sid, lst in ns.items():
            f.write(f"{sid}   {len(lst)}\n")
        f.write("0   # number of side sets\n")
        f.write("# end dimensions\n")

        f.write("*nodesets\n")
        for sid, lst in ns.items():
            f.write("*set\n")
            f.write(f"{len(lst)}   # number of nodes\n")
            for k0 in range(0, len(lst), 10):
                f.write("  ".join(str(n) for n in lst[k0:k0 + 10]) + "\n")
        f.write("# end node sets\n")

        f.write("*sidesets\n")
        f.write("# end side sets\n")

        f.write("*elements\n")
        f.write("*set\n")
        f.write(f"{num_elem}  # number of elements\n")
        f.write("8  # number of element nodes\n")
        for e in elements:
            f.write(
                f"{e[0]}  {e[1]}  {e[2]}  {e[3]}  {e[4]}  "
                f"{e[5]}  {e[6]}  {e[7]}  {e[8]}\n"
            )
        f.write("# end elements\n")

        f.write("*nodes\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("3   # number of spatial dimensions\n")
        for n in nodes:
            f.write(f"    {n[0]}   {n[1]:.10e}   {n[2]:.10e}   {n[3]:.10e}\n")

    print(f"Wrote {out_path}")
    print(f"  Nodes: {num_nodes}, Elements: {num_elem}")
    for sid, lst in ns.items():
        labels = {
            1: "x=0  clamped end",
            2: "x=Lx clamped end",
            3: "y=0  side roller",
            4: "y=Ly side roller",
            5: "z=0  bottom electrode",
            6: "z=Lz top electrode",
        }
        print(f"  NS {sid} ({labels[sid]}): {len(lst)} nodes")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Lx", type=float, default=20.0)
    ap.add_argument("--Ly", type=float, default=2.0)
    ap.add_argument("--Lz", type=float, default=1.0)
    ap.add_argument("--nx", type=int, default=80)
    ap.add_argument("--ny", type=int, default=8)
    ap.add_argument("--nz", type=int, default=4)
    ap.add_argument("--imperfection", type=float, default=0.01,
                    help="amplitude of mode-1 sinusoidal z-offset (default 1%% of Lz)")
    ap.add_argument("--out", default="plate.geom")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out_path = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)
    generate(args.Lx, args.Ly, args.Lz,
             args.nx, args.ny, args.nz,
             args.imperfection, out_path)
