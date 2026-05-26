#!/usr/bin/env python3
"""
Generate a 2D quad4 bar mesh for the crease/wrinkle dielectric-elastomer test,
matching Seifi & Park 2016 IJSS (eq. 4.2 "Surface creasing in constrained 2D
strip"): 80 x 4 DE film with unit-edge square elements.

Domain:  x in [0, Lx],  y in [0, H]
  Defaults Lx = 80, H = 4 (paper geometry), Nx = 80, Ny = 4 (square unit elements)

Mechanics BCs (applied in the XML via node sets):
  bottom (y=0):     clamped, u_x = u_y = 0
  left   (x=0):     roller,  u_x = 0
  right  (x=Lx):    roller,  u_x = 0
  top    (y=H):     free; surface tension acts here

Electrical BCs:
  bottom (y=0):     ground electrode, Psi = 0
  top    (y=H):     driven electrode, Psi = V(t)

Node sets:
  1  y = 0    bottom (ground + mechanical clamp)
  2  y = H    top    (driven electrode; surface-tension side-set lives on faces of these elements)
  3  x = 0    left   (roller)
  4  x = Lx   right  (roller)

Side sets (used to attach surface tension via the
`updated_lagrangian_Q1P0_surface` element's <surface_tension side_set_ID="1" .../>):
  1  top face (y = H) of the elements adjacent to the top edge

Node numbering: i-outer, j-inner.  nid(i, j) = i*(Ny+1) + j + 1
  - i: x-index (0..Nx)
  - j: y-index (0..Ny)

Tahoe quad4 local face IDs (counter-clockwise from bottom):
  face 1 = (n0,n1) bottom   face 2 = (n1,n2) right
  face 3 = (n2,n3) top      face 4 = (n3,n0) left
"""

import argparse
import os


def generate(Lx, H, Nx, Ny, out_path):
    x0 = 0.0
    y0 = 0.0
    nnx, nny = Nx + 1, Ny + 1
    num_nodes = nnx * nny
    num_elem  = Nx * Ny

    def nid(i, j):
        return i * nny + j + 1

    # Nodes — perfectly regular grid, i-outer
    nodes = []
    for i in range(nnx):
        x = x0 + i * Lx / Nx
        for j in range(nny):
            y = y0 + j * H / Ny
            nodes.append((nid(i, j), x, y))

    # Elements (quad4 CCW: bottom-left, bottom-right, top-right, top-left)
    elements = []
    eid_of = {}
    eid = 1
    for j in range(Ny):
        for i in range(Nx):
            conn = (
                nid(i,   j  ), nid(i+1, j  ),
                nid(i+1, j+1), nid(i,   j+1),
            )
            elements.append((eid, *conn))
            eid_of[(i, j)] = eid
            eid += 1

    # Node sets
    ns = {sid: [] for sid in range(1, 5)}
    for i in range(nnx):
        for j in range(nny):
            n = nid(i, j)
            if j == 0:    ns[1].append(n)   # bottom
            if j == Ny:   ns[2].append(n)   # top
            if i == 0:    ns[3].append(n)   # left
            if i == Nx:   ns[4].append(n)   # right
    for sid in ns:
        ns[sid].sort()

    ns_labels = {
        1: "y=0    bottom  (ground electrode + u_x=u_y=0)",
        2: f"y={H:g}  top     (driven electrode; surface tension on face 3)",
        3: "x=0    left    (u_x=0 roller)",
        4: f"x={Lx:g}    right   (u_x=0 roller)",
    }

    # Side sets — only the top edge (face 3 of the elements in the topmost row j=Ny-1)
    side_sets = {
        1: [(eid_of[(i, Ny - 1)], 3) for i in range(Nx)],
    }
    ss_labels = {1: f"top edge (y={H:g}) — surface tension"}

    with open(out_path, "w") as f:
        f.write("*version\n1.0\n")
        f.write("*title\n")
        f.write(f"crease/wrinkle bar: {Nx}x{Ny} quad4, Lx={Lx} H={H}\n")
        f.write("*dimensions\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("2   # number of spatial dimensions\n")
        f.write("1   # number of element sets\n")
        f.write("# [ID] [nel] [nen]\n")
        f.write(f"1  {num_elem}  4\n")
        f.write(f"{len(ns)}   # number of node sets\n")
        f.write("# [ID] [nnd]\n")
        for sid in sorted(ns):
            f.write(f"{sid}  {len(ns[sid])}\n")
        f.write(f"{len(side_sets)}  # number of side sets\n")
        f.write("# [ID] [associated block ID] [n faces]\n")
        for sid in sorted(side_sets):
            f.write(f"{sid}  1  {len(side_sets[sid])}\n")
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

        f.write("*elements\n*set\n")
        f.write(f"{num_elem}  # number of elements\n")
        f.write("4  # number of element nodes\n")
        for e in elements:
            f.write("  " + "  ".join(str(v) for v in e) + "\n")
        f.write("# end elements\n")

        f.write("*nodes\n")
        f.write(f"{num_nodes}  # number of nodes\n")
        f.write("2   # number of spatial dimensions\n")
        for nd in nodes:
            f.write(f"  {nd[0]}  {nd[1]:.10e}  {nd[2]:.10e}\n")

    print(f"Written: {out_path}")
    print(f"  Nodes: {num_nodes}   Elements: {num_elem}")
    for sid in sorted(ns):
        print(f"  NS {sid} ({ns_labels[sid]}): {len(ns[sid])} nodes")
    for sid in sorted(side_sets):
        print(f"  SS {sid} ({ss_labels[sid]}): {len(side_sets[sid])} faces")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lx", type=float, default=80.0,
                    help="bar length (default 80, matches Seifi-Park 2016 IJSS)")
    ap.add_argument("--H",  type=float, default=4.0,
                    help="bar thickness (default 4, matches paper)")
    ap.add_argument("--Nx", type=int,   default=80,
                    help="elements along x (default 80 -> dx = Lx/Nx = 1)")
    ap.add_argument("--Ny", type=int,   default=4,
                    help="elements through y-thickness (default 4 -> dy = H/Ny = 1)")
    ap.add_argument("--out", default="../meshes/bar_2D.geom")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)
    generate(args.Lx, args.H, args.Nx, args.Ny, out)
