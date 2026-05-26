#!/usr/bin/env python3
"""
Generate a 3D hex8 bar mesh for the crease/wrinkle dielectric-elastomer test,
3D extension of the 2D Seifi & Park 2016 IJSS replication (section 4.2).

Domain:  x in [0, Lx],  y in [0, H],  z in [0, Lz]
  Defaults Lx = 80, H = 4, Lz = 4  (unit-cube hex elements)
           Nx = 80, Ny = 4, Nz = 4

With z-direction rollers on the front (z=0) and back (z=Lz) faces, the 3D
problem reduces to plane strain in the x-y plane and reproduces the 2D
crease/wrinkle results.  Relax those z-rollers to study true 3D mode
selection (2D wrinkle wavevector kx, kz in the top-face plane).

Mechanics BCs (applied in the XML via node sets):
  bottom (y=0):     clamped, u_x = u_y = u_z = 0
  left   (x=0):     roller,  u_x = 0
  right  (x=Lx):    roller,  u_x = 0
  front  (z=0):     roller,  u_z = 0  (plane-strain mode)
  back   (z=Lz):    roller,  u_z = 0  (plane-strain mode)
  top    (y=H):     free; surface tension acts here

Electrical BCs:
  bottom (y=0):     ground electrode, Psi = 0
  top    (y=H):     driven electrode, Psi = V(t)

Node sets:
  1  y = 0   bottom (ground + mechanical clamp)
  2  y = H   top    (driven electrode; surface tension via side set 1)
  3  x = 0   left   (D_X = 0 roller)
  4  x = Lx  right  (D_X = 0 roller)
  5  z = 0   front  (D_Z = 0 roller)
  6  z = Lz  back   (D_Z = 0 roller)

Side sets:
  1  top face (y = H) of every element in the top row -> 80*Nz = 320 faces
     for the default 80x4x4 mesh.  Used by the surface_tension sub-list of
     the updated_lagrangian_Q1P0_3D_surface element.

Node numbering: i-outer, j-mid, k-inner.
  nid(i, j, k) = i * (Ny+1) * (Nz+1) + j * (Nz+1) + k + 1
  - i: x-index (0..Nx)
  - j: y-index (0..Ny)
  - k: z-index (0..Nz)

Tahoe hex8 local face IDs (1-indexed in the geom file).  From
HexahedronT::NodesOnFacet the (0-indexed) facets are:
  facet 0 (geom 1)  z=z_min  bottom    nodes 0,3,2,1
  facet 1 (geom 2)  z=z_max  top z+    nodes 4,5,6,7
  facet 2 (geom 3)  y=y_min  front y-  nodes 0,1,5,4
  facet 3 (geom 4)  x=x_max  right x+  nodes 1,2,6,5
  facet 4 (geom 5)  y=y_max  TOP y+    nodes 2,3,7,6   <-- surface tension here
  facet 5 (geom 6)  x=x_min  left  x-  nodes 3,0,4,7
"""

import argparse
import os


def generate(Lx, H, Lz, Nx, Ny, Nz, out_path):
    nnx, nny, nnz = Nx + 1, Ny + 1, Nz + 1
    num_nodes = nnx * nny * nnz
    num_elem  = Nx * Ny * Nz

    def nid(i, j, k):
        return i * nny * nnz + j * nnz + k + 1

    # Nodes — perfectly regular grid, i-outer
    nodes = []
    for i in range(nnx):
        x = i * Lx / Nx
        for j in range(nny):
            y = j * H / Ny
            for k in range(nnz):
                z = k * Lz / Nz
                nodes.append((nid(i, j, k), x, y, z))

    # Elements: hex8 with local node order matching Tahoe's HexahedronT.
    # Local node 0..7 placed at:
    #   0 = (i,   j,   k  )    4 = (i,   j,   k+1)
    #   1 = (i+1, j,   k  )    5 = (i+1, j,   k+1)
    #   2 = (i+1, j+1, k  )    6 = (i+1, j+1, k+1)
    #   3 = (i,   j+1, k  )    7 = (i,   j+1, k+1)
    # Facet 4 (y+) = nodes 2,3,7,6 — all at y = (j+1)*H/Ny, i.e. y=H for j=Ny-1.
    elements = []
    eid_of = {}
    eid = 1
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                conn = (
                    nid(i,   j,   k  ),
                    nid(i+1, j,   k  ),
                    nid(i+1, j+1, k  ),
                    nid(i,   j+1, k  ),
                    nid(i,   j,   k+1),
                    nid(i+1, j,   k+1),
                    nid(i+1, j+1, k+1),
                    nid(i,   j+1, k+1),
                )
                elements.append((eid, *conn))
                eid_of[(i, j, k)] = eid
                eid += 1

    # Node sets
    ns = {sid: [] for sid in range(1, 7)}
    for i in range(nnx):
        for j in range(nny):
            for k in range(nnz):
                n = nid(i, j, k)
                if j == 0:    ns[1].append(n)  # bottom
                if j == Ny:   ns[2].append(n)  # top
                if i == 0:    ns[3].append(n)  # left
                if i == Nx:   ns[4].append(n)  # right
                if k == 0:    ns[5].append(n)  # front
                if k == Nz:   ns[6].append(n)  # back
    for sid in ns:
        ns[sid].sort()

    ns_labels = {
        1: "y=0    bottom  (ground electrode + u_x=u_y=u_z=0)",
        2: f"y={H:g}  top     (driven electrode; surface tension on face 5)",
        3: "x=0    left    (u_x=0 roller)",
        4: f"x={Lx:g}   right   (u_x=0 roller)",
        5: "z=0    front   (u_z=0 roller; plane-strain mode)",
        6: f"z={Lz:g}    back    (u_z=0 roller; plane-strain mode)",
    }

    # Side sets — top face (y=H) of elements in the top row j=Ny-1
    # Local face id for y+ is facet 4 (0-indexed) -> geom face 5 (1-indexed).
    side_sets = {
        1: [(eid_of[(i, Ny - 1, k)], 5) for i in range(Nx) for k in range(Nz)],
    }
    ss_labels = {1: f"top face (y={H:g}) — surface tension"}

    with open(out_path, "w") as f:
        f.write("*version\n1.0\n")
        f.write("*title\n")
        f.write(f"crease/wrinkle bar 3D: {Nx}x{Ny}x{Nz} hex8, Lx={Lx} H={H} Lz={Lz}\n")
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
    for sid in sorted(side_sets):
        print(f"  SS {sid} ({ss_labels[sid]}): {len(side_sets[sid])} faces")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Lx", type=float, default=80.0,
                    help="bar length (default 80, matches Seifi-Park 2016)")
    ap.add_argument("--H",  type=float, default=4.0,
                    help="bar thickness (default 4, matches paper)")
    ap.add_argument("--Lz", type=float, default=4.0,
                    help="bar width in z (default 4 -> unit cube cells)")
    ap.add_argument("--Nx", type=int,   default=80,
                    help="elements along x (default 80 -> dx = Lx/Nx = 1)")
    ap.add_argument("--Ny", type=int,   default=4,
                    help="elements through y-thickness (default 4 -> dy = 1)")
    ap.add_argument("--Nz", type=int,   default=4,
                    help="elements along z (default 4 -> dz = 1)")
    ap.add_argument("--out", default="../meshes/bar_3D.geom")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)
    generate(args.Lx, args.H, args.Lz, args.Nx, args.Ny, args.Nz, out)
