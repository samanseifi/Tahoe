#!/usr/bin/env python3
"""Generate a 2D quad4 plate mesh for dielectric-elastomer buckling (plane-strain).

Domain (default): x in [0,Lx], y in [0,Ly]
  Lx = 20  long axis (slender direction)
  Ly =  1  thickness  -- electric field applied across this dimension

The plate acts as a beam compressed by the in-plane Maxwell stress
  sigma_E_xx = -0.5 * epsilon * (V/Ly)^2
which drives Euler buckling in the transverse (y) direction when sigma_E_xx
exceeds the linearized critical load sigma_cr ~ 4*pi^2*mu/(3*(Ly/Lx)^2).

A small mode-1 sinusoidal Y imperfection seeds the instability:
  y_new = y + amp * sin(pi * x / Lx)   (vanishes at clamped ends)

Node sets:
  1  x = 0     left  clamped end (fix D_X, D_Y)
  2  x = Lx    right clamped end (fix D_X, D_Y)
  3  y = 0     bottom electrode (Psi = 0)
  4  y = Ly    top    electrode (Psi = V(t))

Output:
  plate_2D.geom  -- Tahoe single-file ASCII geometry
"""
import argparse
import math
import os


def generate(Lx, Ly, nx, ny, imperfection, out_path):
    nnx, nny = nx + 1, ny + 1
    num_nodes = nnx * nny
    num_elem = nx * ny

    def nid(i, j):
        # row-major: row j (bottom=0), column i (left=0) → 1-indexed
        return j * nnx + i + 1

    nodes = []
    for j in range(nny):
        for i in range(nnx):
            x = i * Lx / nx
            y = j * Ly / ny
            y += imperfection * math.sin(math.pi * x / Lx)
            nodes.append((nid(i, j), x, y))

    elements = []
    eid = 1
    for j in range(ny):
        for i in range(nx):
            bl = nid(i,     j)
            br = nid(i + 1, j)
            tr = nid(i + 1, j + 1)
            tl = nid(i,     j + 1)
            elements.append((eid, bl, br, tr, tl))
            eid += 1

    ns = {1: [], 2: [], 3: [], 4: []}
    for j in range(nny):
        for i in range(nnx):
            n = nid(i, j)
            if i == 0:  ns[1].append(n)
            if i == nx: ns[2].append(n)
            if j == 0:  ns[3].append(n)
            if j == ny: ns[4].append(n)

    with open(out_path, "w") as f:
        f.write("*version\n1.0\n")
        f.write("*title\n")
        f.write(
            f"Dielectric elastomer plate 2D: {nx}x{ny} quad4, "
            f"Lx={Lx} Ly={Ly}, imperfection={imperfection}\n"
        )

        f.write("*dimensions\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("2   # number of spatial dimensions\n")
        f.write("1   # number of element sets\n")
        f.write("# [ID] [nel] [nen]\n")
        f.write(f"1   {num_elem}   4\n")
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
        f.write("4  # number of element nodes\n")
        for e in elements:
            f.write(f"{e[0]}  {e[1]}  {e[2]}  {e[3]}  {e[4]}\n")
        f.write("# end elements\n")

        f.write("*nodes\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("2   # number of spatial dimensions\n")
        for n in nodes:
            f.write(f"    {n[0]}   {n[1]:.10e}   {n[2]:.10e}\n")

    print(f"Wrote {out_path}")
    print(f"  Nodes: {num_nodes}, Elements: {num_elem}")
    for sid, lst in ns.items():
        labels = {
            1: "x=0  left clamped end",
            2: "x=Lx right clamped end",
            3: "y=0  bottom electrode (Psi=0)",
            4: "y=Ly top electrode    (Psi=V(t))",
        }
        print(f"  NS {sid} ({labels[sid]}): {len(lst)} nodes")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Lx", type=float, default=20.0)
    ap.add_argument("--Ly", type=float, default=1.0)
    ap.add_argument("--nx", type=int, default=80)
    ap.add_argument("--ny", type=int, default=4)
    ap.add_argument("--imperfection", type=float, default=0.01,
                    help="amplitude of mode-1 sinusoidal y-offset (default 1%% of Ly)")
    ap.add_argument("--out", default="plate_2D.geom")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out_path = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)
    generate(args.Lx, args.Ly, args.nx, args.ny, args.imperfection, out_path)
