#!/usr/bin/env python3
"""Generate a Tahoe .geom (ModelManager text) for a Lx x H rectangular film,
Nx x Ny quad4 elements.  Node sets: 1=bottom, 2=top, 3=left, 4=right.
Side set 1 = top surface (for the surface_tension element).  Matches the layout
of film_L80.geom exactly (column-major nodes, row-major elements, top=face 3).
Usage: python3 generate_film.py <Lx> <H> <Nx> <Ny> <out.geom>
"""
import sys

Lx, H = float(sys.argv[1]), float(sys.argv[2])
Nx, Ny = int(sys.argv[3]), int(sys.argv[4])
out = sys.argv[5]
dx, dy = Lx / Nx, H / Ny
nnp = (Nx + 1) * (Ny + 1)
nel = Nx * Ny


def nid(i, j):                       # column-major node id (j fastest)
    return i * (Ny + 1) + j + 1


L = []
L.append("*version\n1.0")
L.append(f"*title\ncrease/wrinkle bar: {Nx}x{Ny} quad4, Lx={Lx} H={H}")
L.append("*dimensions")
L.append(f"{nnp}  # number of nodes")
L.append("2   # number of spatial dimensions")
L.append("1   # number of element sets")
L.append("# [ID] [nel] [nen]")
L.append(f"1  {nel}  4")
L.append("5   # number of node sets")
L.append("# [ID] [nnd]")
L.append(f"1  {Nx+1}")
L.append(f"2  {Nx+1}")
L.append(f"3  {Ny+1}")
L.append(f"4  {Ny+1}")
L.append("5  1")            # NS5 = bottom-left corner (single node, x-pin for periodic BC)
L.append("1  # number of side sets")
L.append("# [ID] [associated block ID] [n faces]")
L.append(f"1  1  {Nx}")
L.append("# end dimensions")

# node sets
def nodeset(ids):
    s = [f"{len(ids)}  # number of nodes"]
    for k in range(0, len(ids), 10):
        s.append("  ".join(str(v) for v in ids[k:k + 10]))
    return "\n".join(s)

bottom = [nid(i, 0) for i in range(Nx + 1)]
top = [nid(i, Ny) for i in range(Nx + 1)]
left = [nid(0, j) for j in range(Ny + 1)]
right = [nid(Nx, j) for j in range(Ny + 1)]
corner = [nid(0, 0)]                     # bottom-left corner
L.append("*nodesets")
for ids in (bottom, top, left, right, corner):
    L.append("*set")
    L.append(nodeset(ids))
L.append("# end node sets")

# side set: top faces (element (i, Ny-1), local face 3)
L.append("*sidesets")
L.append("*set")
L.append(f"{Nx}  # number of faces")
for i in range(Nx):
    eid = (Ny - 1) * Nx + i + 1
    L.append(f"  {eid}  3")
L.append("# end side sets")

# elements (row-major: j outer, i inner; id = j*Nx + i + 1)
L.append("*elements")
L.append("*set")
L.append(f"{nel}  # number of elements")
L.append("4  # number of element nodes")
for j in range(Ny):
    for i in range(Nx):
        eid = j * Nx + i + 1
        n1, n2 = nid(i, j), nid(i + 1, j)
        n3, n4 = nid(i + 1, j + 1), nid(i, j + 1)
        L.append(f"  {eid}  {n1}  {n2}  {n3}  {n4}")
L.append("# end elements")

# nodes (column-major)
L.append("*nodes")
L.append(f"{nnp}  # number of nodes")
L.append("2   # number of spatial dimensions")
for i in range(Nx + 1):
    for j in range(Ny + 1):
        L.append(f"  {nid(i,j)}  {i*dx:.10e}  {j*dy:.10e}")

open(out, "w").write("\n".join(L) + "\n")
print(f"wrote {out}: {nnp} nodes, {nel} elems, dx={dx:.3f} dy={dy:.3f}, "
      f"NS bottom/top={Nx+1} left/right={Ny+1}, SS1 top={Nx} faces")
