#!/usr/bin/env python3
"""
Generate a 2D quad4 *bilayer* bar mesh for the multilayer crease/wrinkle
dielectric-elastomer study (interfacial surface energy extension of
Seifi & Park 2016 IJSS).

Two bonded layers stacked in y, sharing a conforming row of nodes on the
internal interface:

    y in [H1, H1+H2]   block 2  (top layer / "film")    shear modulus mu2
    --------------------  internal interface y = H1  (gamma_int acts here)
    y in [0,  H1   ]   block 1  (bottom layer / "substrate") shear modulus mu1
    ====================  rigid base y = 0

The two blocks share the interface nodes (displacement-continuous bond),
so a side set placed on the top faces (face 3) of the topmost substrate
elements coincides geometrically with the bottom faces of the film and
applies the interfacial Young-Laplace / Gurtin-Murdoch traction there.
SimoQ1P0_Surface merges these interior-interface faces into its surface
list (see SimoQ1P0_Surface.cpp ~L162-216).

Node sets:
  1  y = 0        bottom (ground electrode + mechanical clamp u_x=u_y=0)
  2  y = H1+H2    top    (driven electrode; top surface tension)
  3  x = 0        left   (roller u_x=0)
  4  x = Lx       right  (roller u_x=0)
  5  y = H1       interface nodes (diagnostics / optional BC)

Side sets:
  1  top face  (y = H1+H2)  -> top surface tension gamma_top
  2  interface (y = H1)     -> interfacial surface tension gamma_int

Block IDs:
  1  substrate (j = 0 .. Ny1-1)
  2  film      (j = Ny1 .. Ny1+Ny2-1)

Node numbering (single global grid, i-outer): nid(i,j) = i*nny + j + 1
with nny = Ny1 + Ny2 + 1 rows of nodes.
"""

import argparse
import os


def generate(Lx, H1, H2, Nx, Ny1, Ny2, out_path):
    nny = Ny1 + Ny2 + 1
    nnx = Nx + 1
    num_nodes = nnx * nny
    Ny = Ny1 + Ny2
    num_elem = Nx * Ny

    def nid(i, j):
        return i * nny + j + 1

    # y-coordinate of node row j: layer 1 occupies rows 0..Ny1, layer 2 Ny1..Ny
    def yj(j):
        if j <= Ny1:
            return H1 * j / Ny1
        return H1 + H2 * (j - Ny1) / Ny2

    # Nodes
    nodes = []
    for i in range(nnx):
        x = i * Lx / Nx
        for j in range(nny):
            nodes.append((nid(i, j), x, yj(j)))

    # Elements, split into two blocks by j.  Tahoe .geom convention: element
    # IDs restart at 1 (LOCAL) within each block; side sets reference block-local
    # element IDs.  Track per-block local id for each (i,j).
    elems_b1, elems_b2 = [], []
    lid_of = {}   # (i,j) -> (block, local_id)
    lid1 = lid2 = 0
    for j in range(Ny):
        for i in range(Nx):
            conn = (nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
            if j < Ny1:
                lid1 += 1
                elems_b1.append((lid1, *conn))
                lid_of[(i, j)] = (1, lid1)
            else:
                lid2 += 1
                elems_b2.append((lid2, *conn))
                lid_of[(i, j)] = (2, lid2)

    # Node sets
    ns = {sid: [] for sid in range(1, 6)}
    for i in range(nnx):
        for j in range(nny):
            n = nid(i, j)
            if j == 0:    ns[1].append(n)   # bottom
            if j == Ny:   ns[2].append(n)   # top
            if i == 0:    ns[3].append(n)   # left
            if i == Nx:   ns[4].append(n)   # right
            if j == Ny1:  ns[5].append(n)   # interface
    for sid in ns:
        ns[sid].sort()

    # Side sets reference block-LOCAL element IDs.
    #   top surface = face 3 of top film row (j=Ny-1)      -> block 2
    #   interface   = face 3 of top substrate row (j=Ny1-1) -> block 1
    side_sets = {
        1: [(lid_of[(i, Ny - 1)][1], 3) for i in range(Nx)],
        2: [(lid_of[(i, Ny1 - 1)][1], 3) for i in range(Nx)],
    }
    # side set 2 lives on block 1 elements; side set 1 on block 2 elements
    ss_block = {1: 2, 2: 1}

    with open(out_path, "w") as f:
        f.write("*version\n1.0\n")
        f.write("*title\n")
        f.write(f"bilayer: {Nx}x{Ny} quad4, Lx={Lx} H1={H1} H2={H2}\n")
        f.write("*dimensions\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("2   # number of spatial dimensions\n")
        f.write("2   # number of element sets\n")
        f.write("# [ID] [nel] [nen]\n")
        f.write(f"1  {len(elems_b1)}  4\n")
        f.write(f"2  {len(elems_b2)}  4\n")
        f.write(f"{len(ns)}   # number of node sets\n")
        f.write("# [ID] [nnd]\n")
        for sid in sorted(ns):
            f.write(f"{sid}  {len(ns[sid])}\n")
        f.write(f"{len(side_sets)}  # number of side sets\n")
        f.write("# [ID] [associated block ID] [n faces]\n")
        for sid in sorted(side_sets):
            f.write(f"{sid}  {ss_block[sid]}  {len(side_sets[sid])}\n")
        f.write("# end dimensions\n")

        f.write("*nodesets\n")
        for sid in sorted(ns):
            lst = ns[sid]
            f.write("*set\n")
            f.write(f"{len(lst)}  # number of nodes\n")
            for k0 in range(0, len(lst), 10):
                f.write("  ".join(str(n) for n in lst[k0:k0 + 10]) + "\n")
        f.write("# end node sets\n")

        f.write("*sidesets\n")
        for sid in sorted(side_sets):
            faces = side_sets[sid]
            f.write("*set\n")
            f.write(f"{len(faces)}  # number of faces\n")
            for (e, fc) in faces:
                f.write(f"  {e}  {fc}\n")
        f.write("# end side sets\n")

        f.write("*elements\n")
        for blk_id, elems in ((1, elems_b1), (2, elems_b2)):
            f.write("*set\n")
            f.write(f"{len(elems)}  # number of elements\n")
            f.write("4  # number of element nodes\n")
            for e in elems:
                f.write("  " + "  ".join(str(v) for v in e) + "\n")
        f.write("# end elements\n")

        f.write("*nodes\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("2   # number of spatial dimensions\n")
        for nd in nodes:
            f.write(f"  {nd[0]}  {nd[1]:.10e}  {nd[2]:.10e}\n")

    print(f"Written: {out_path}")
    print(f"  Nodes: {num_nodes}  Elements: {num_elem} (block1={len(elems_b1)}, block2={len(elems_b2)})")
    print(f"  Interface at y={H1:g} ({len(ns[5])} nodes), total H={H1+H2:g}")
    for sid in sorted(ns):
        print(f"  NS {sid}: {len(ns[sid])} nodes")
    for sid in sorted(side_sets):
        print(f"  SS {sid} (block {ss_block[sid]}): {len(side_sets[sid])} faces")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lx", type=float, default=80.0)
    ap.add_argument("--H1", type=float, default=2.0, help="substrate thickness")
    ap.add_argument("--H2", type=float, default=2.0, help="film thickness")
    ap.add_argument("--Nx", type=int, default=160, help="elements along x")
    ap.add_argument("--Ny1", type=int, default=8, help="elements through substrate")
    ap.add_argument("--Ny2", type=int, default=8, help="elements through film")
    ap.add_argument("--out", default="../meshes/bilayer_2D.geom")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)
    generate(args.Lx, args.H1, args.H2, args.Nx, args.Ny1, args.Ny2, out)
