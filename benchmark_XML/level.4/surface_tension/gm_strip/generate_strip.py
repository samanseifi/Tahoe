#!/usr/bin/env python3
"""
Tiny 8x2 quad4 strip for the #54 Gurtin-Murdoch benchmark.

Domain: x in [0, 8], y in [0, 1].  Plane strain.

Node sets:
  1  y = 0    bottom (D_Y = 0, substrate-bonded)
  2  y = 1    top    (free surface; surface_tension lives on side set 1)
  3  x = 0    left   (D_X = 0)
  4  x = 8    right  (D_X = prescribed stretch)

Side set:
  1  top edge (y = 1), 8 face-3 (Tahoe 1-indexed) entries
"""
import os

Lx, H = 8.0, 1.0
Nx, Ny = 8, 2
nnx, nny = Nx + 1, Ny + 1

def nid(i, j):
    return i * nny + j + 1   # i-outer, j-inner

nodes = []
for i in range(nnx):
    for j in range(nny):
        nodes.append((nid(i, j), i * Lx / Nx, j * H / Ny))

# elements, CCW quad4
elements = []
eid_of = {}
eid = 1
for j in range(Ny):
    for i in range(Nx):
        elements.append((eid,
                         nid(i,   j  ),
                         nid(i+1, j  ),
                         nid(i+1, j+1),
                         nid(i,   j+1)))
        eid_of[(i, j)] = eid
        eid += 1

ns = {1: [], 2: [], 3: [], 4: []}
for i in range(nnx):
    for j in range(nny):
        n = nid(i, j)
        if j == 0:    ns[1].append(n)
        if j == Ny:   ns[2].append(n)
        if i == 0:    ns[3].append(n)
        if i == Nx:   ns[4].append(n)
for k in ns: ns[k].sort()

side_sets = {1: [(eid_of[(i, Ny-1)], 3) for i in range(Nx)]}

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "strip.geom")
with open(out, "w") as f:
    f.write("*version\n1.0\n*title\nGM strip benchmark %dx%d\n" % (Nx, Ny))
    f.write("*dimensions\n")
    f.write(f"{len(nodes)}  # number of nodes\n2   # spatial dim\n1   # element sets\n# [ID] [nel] [nen]\n1  {len(elements)}  4\n")
    f.write(f"{len(ns)}   # node sets\n# [ID] [nnd]\n")
    for k in sorted(ns): f.write(f"{k}  {len(ns[k])}\n")
    f.write(f"{len(side_sets)}  # side sets\n# [ID] [block] [n_faces]\n")
    for k in sorted(side_sets): f.write(f"{k}  1  {len(side_sets[k])}\n")
    f.write("# end dimensions\n")
    f.write("*nodesets\n")
    for k in sorted(ns):
        f.write(f"*set\n{len(ns[k])}\n")
        for n in ns[k]: f.write(f"  {n}\n")
    f.write("# end node sets\n")
    f.write("*sidesets\n")
    for k in sorted(side_sets):
        f.write(f"*set\n{len(side_sets[k])}\n")
        for e, fc in side_sets[k]: f.write(f"  {e}  {fc}\n")
    f.write("# end side sets\n")
    f.write("*elements\n*set\n")
    f.write(f"{len(elements)}\n4\n")
    for e in elements:
        f.write("  " + "  ".join(str(v) for v in e) + "\n")
    f.write("# end elements\n*nodes\n")
    f.write(f"{len(nodes)}\n2\n")
    for n in nodes:
        f.write(f"  {n[0]}  {n[1]:.10e}  {n[2]:.10e}\n")

print(f"Written: {out}  ({len(nodes)} nodes, {len(elements)} elements, "
      f"top side set: {len(side_sets[1])} faces)")
