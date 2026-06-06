#!/usr/bin/env python3
"""Generate a structured hex cube for the meshfree 90-degree rotation example.

Cube:  x in [0, L]  (x=0 face fixed, x=L face rotated),
       y,z in [-L/2, +L/2]  (centered so the rotation axis is the x-axis).

Emits Tahoe geometry: cube.geom (index) + cube.node + cube.elem.
Node ordering: i (x) fastest, then j (y), then k (z).
Node sets: 1 = fixed face (x=0), 2 = rotating face (x=L).
"""
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4   # elements per edge
L = 1.0
base = "cube"

npe = N + 1
def nid(i, j, k):                # 1-based node id
    return i + j * npe + k * npe * npe + 1

# ---- nodes ----
coords = []
for k in range(npe):
    for j in range(npe):
        for i in range(npe):
            x = i * L / N
            y = -L / 2 + j * L / N
            z = -L / 2 + k * L / N
            coords.append((nid(i, j, k), x, y, z))
coords.sort()

with open(f"{base}.node", "w") as f:
    f.write(f"{len(coords)}  # number of nodes\n3   # number of spatial dimensions\n")
    for (n, x, y, z) in coords:
        f.write(f"{n:6d}  {x: .7e}  {y: .7e}  {z: .7e}\n")

# ---- elements (hex8, Tahoe ordering: bottom CCW then top CCW) ----
elems = []
eid = 0
for k in range(N):
    for j in range(N):
        for i in range(N):
            eid += 1
            n = [nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                 nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1)]
            elems.append((eid, n))

with open(f"{base}.elem", "w") as f:
    f.write(f"{len(elems)}  # number of elements\n8   # number of element nodes\n")
    for (e, n) in elems:
        f.write(f"{e:6d}  " + "  ".join(f"{m}" for m in n) + "\n")

# ---- node sets ----
fixed   = [nid(0, j, k) for k in range(npe) for j in range(npe)]   # x = 0
rotate  = [nid(N, j, k) for k in range(npe) for j in range(npe)]   # x = L

def emit_set(f, ids):
    f.write(f"{len(ids)}   # number of nodes\n")
    for c in range(0, len(ids), 8):
        f.write(" ".join(f"{m:4d}" for m in ids[c:c+8]) + "\n")

# ---- geom index ----
with open(f"{base}.geom", "w") as f:
    f.write("*version\n1.0\n*title\n")
    f.write(f"{N}x{N}x{N} hex cube, centered in y,z; x=0 fixed, x=L rotated\n")
    f.write("*dimensions\n")
    f.write(f"{len(coords)}  # number of nodes\n3   # number of spatial dimensions\n")
    f.write("1   # number of element sets\n# [ID] [nel] [nen]\n")
    f.write(f"1   {len(elems)}   8\n")
    f.write("2   # number of node sets\n# [ID] [nnd]\n")
    f.write(f"1   {len(fixed)}\n2   {len(rotate)}\n")
    f.write("0   # number of side sets\n# end dimensions\n")
    f.write("*nodesets\n")
    f.write("*set\n");  emit_set(f, fixed)
    f.write("*set\n");  emit_set(f, rotate)
    f.write("# end node sets\n")
    f.write("*sidesets\n")
    f.write("*elements\n*set\n")
    f.write(f"{base}.elem\n# end elements\n")
    f.write("*nodes\n")
    f.write(f"{base}.node\n")

print(f"generated {base}.geom/.node/.elem : N={N}, {len(coords)} nodes, {len(elems)} elements")
print(f"fixed-face nodes (set 1): {len(fixed)}; rotating-face nodes (set 2): {len(rotate)}")
