#!/usr/bin/env python3
"""Generate a two-stacked-cubes mesh with refined Hex8 grids and contact
sidesets — INLINE Tahoe geom format (single .geom file, no external
.ns / .ss / .es / .nd files).  Inline is what cube.0.geom uses and is
the format the legacy parser handles for multi-block + sidesets.

Two cubes:
  cube 1 (block 1):  (x,y,z) in [-0.5, 0.5]^2 x [0, 1]
  cube 2 (block 2):  same in xy, z in [1, 2]
Each cube discretised as N x N x N Hex8 cells.

Node-set IDs:
  1 — bottom face of cube 1 (z=0): clamped
  2 — top face of cube 2 (z=2): driven
  3 — interface nodes belonging to cube 1 (z=1, cube-1 side)
  4 — interface nodes belonging to cube 2 (z=1, cube-2 side)

Side-set IDs (for contact):
  1 — top face of cube 1 (block 1, face 2 = +z)
  2 — bottom face of cube 2 (block 2, face 0 = -z)

Element block IDs:
  1 — cube 1 elements (block 1 local IDs 1..n^3)
  2 — cube 2 elements (block 2 local IDs 1..n^3)
"""
import os


def write_mesh(name: str, n: int = 5):
    nn = n + 1
    cube1_nnod = nn * nn * nn
    cube2_nnod = nn * nn * nn
    numnod = cube1_nnod + cube2_nnod
    cube_nelem = n * n * n
    numel = 2 * cube_nelem

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)

    def nid_c1(i, j, k):
        return k * nn * nn + j * nn + i + 1
    def nid_c2(i, j, k):
        return cube1_nnod + k * nn * nn + j * nn + i + 1

    # Node coords
    coords = []
    for k in range(nn):
        for j in range(nn):
            for i in range(nn):
                coords.append((nid_c1(i, j, k),
                               -0.5 + i / n, -0.5 + j / n, 0.0 + k / n))
    for k in range(nn):
        for j in range(nn):
            for i in range(nn):
                coords.append((nid_c2(i, j, k),
                               -0.5 + i / n, -0.5 + j / n, 1.0 + k / n))

    # Element connectivities (block 1 then block 2)
    conn_c1 = []
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                ns = [nid_c1(i,   j,   k),   nid_c1(i+1, j,   k),
                      nid_c1(i+1, j+1, k),   nid_c1(i,   j+1, k),
                      nid_c1(i,   j,   k+1), nid_c1(i+1, j,   k+1),
                      nid_c1(i+1, j+1, k+1), nid_c1(i,   j+1, k+1)]
                conn_c1.append((eid, ns)); eid += 1

    conn_c2 = []
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                ns = [nid_c2(i,   j,   k),   nid_c2(i+1, j,   k),
                      nid_c2(i+1, j+1, k),   nid_c2(i,   j+1, k),
                      nid_c2(i,   j,   k+1), nid_c2(i+1, j,   k+1),
                      nid_c2(i+1, j+1, k+1), nid_c2(i,   j+1, k+1)]
                conn_c2.append((eid, ns)); eid += 1

    # Node sets
    ns_bottom = [nid_c1(i, j, 0) for j in range(nn) for i in range(nn)]
    ns_top    = [nid_c2(i, j, n) for j in range(nn) for i in range(nn)]
    ns_iface1 = [nid_c1(i, j, n) for j in range(nn) for i in range(nn)]
    ns_iface2 = [nid_c2(i, j, 0) for j in range(nn) for i in range(nn)]

    # Side sets — element IDs LOCAL to their associated block, face IDs are
    # 1-INDEXED in Tahoe geom format:
    #   face 1 = bottom (-z, dat8 row 0: nodes 0,3,2,1)
    #   face 2 = top    (+z, dat8 row 1: nodes 4,5,6,7)
    ss1 = []  # cube 1 top layer (k=n-1), face 2 (+z)
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                if k == n - 1: ss1.append((eid, 2))
                eid += 1
    ss2 = []  # cube 2 bottom layer (k=0), face 1 (-z)
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                if k == 0: ss2.append((eid, 1))
                eid += 1

    # Write inline geom file
    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\ntwo_stacked_cubes_refined\n")
        f.write("*dimensions\n")
        f.write(f"{numnod}   # number of nodes\n")
        f.write("3   # number of spatial dimensions\n")
        f.write("2   # number of element sets\n")
        f.write("# [ID] [nel] [nen]\n")
        f.write(f"1   {cube_nelem}   8\n")
        f.write(f"2   {cube_nelem}   8\n")
        f.write("4   # number of node sets\n")
        f.write("# [ID] [nnd]\n")
        f.write(f"1   {len(ns_bottom)}\n")
        f.write(f"2   {len(ns_top)}\n")
        f.write(f"3   {len(ns_iface1)}\n")
        f.write(f"4   {len(ns_iface2)}\n")
        f.write("2   # number of side sets\n")
        f.write("# [ID] [assoc elem block] [ns]\n")
        f.write(f"1  1  {len(ss1)}\n")
        f.write(f"2  2  {len(ss2)}\n")
        f.write("# end dimensions\n")

        f.write("*nodesets\n")
        for lst in [ns_bottom, ns_top, ns_iface1, ns_iface2]:
            f.write("*set\n")
            f.write(f"{len(lst)}   # number of nodes\n")
            f.write(" ".join(str(n) for n in lst) + "\n")
        f.write("# end node sets\n")

        f.write("*sidesets\n")
        for lst in [ss1, ss2]:
            f.write("*set\n")
            f.write(f"{len(lst)}\n")
            for elem, face in lst:
                f.write(f"{elem}  {face}\n")
        f.write("# end side sets\n")

        f.write("*elements\n")
        for conn in [conn_c1, conn_c2]:
            f.write("*set\n")
            f.write(f"{len(conn)}   # number of elements\n")
            f.write("8   # number of element nodes\n")
            for eid, ns in conn:
                f.write(f"{eid}  " + "  ".join(str(x) for x in ns) + "\n")
        f.write("# end elements\n")

        f.write("*nodes\n")
        f.write(f"{numnod}   # number of nodes\n")
        f.write("3   # number of spatial dimensions\n")
        for nid, x, y, z in coords:
            f.write(f"{nid}  {x:.15e}  {y:.15e}  {z:.15e}\n")

    # Clean up any stale external files
    for ext in ("nd", "es0", "es1", "ns0", "ns1", "ns2", "ns3", "ss0", "ss1"):
        p = f"{base}.geom.{ext}"
        if os.path.exists(p):
            os.remove(p)

    print(f"{name}: {n}x{n}x{n} per cube => {numel} hex, {numnod} nodes, INLINE format")


if __name__ == "__main__":
    write_mesh("two_cubes_refined", n=5)
