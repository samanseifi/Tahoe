#!/usr/bin/env python3
"""Generate a mesh for coupled mechanics + phase-field test.

A rectangular bar in 2D (plane strain) under uniaxial tension.
Phase-field d=1 is prescribed at the center to simulate a crack.

Domain: x in [0, 2], y in [0, 1]
Mesh: 40 x 20 quad4 elements
"""

import numpy as np

Lx = 2.0
Ly = 1.0
nx = 40
ny = 20

x_coords = np.linspace(0, Lx, nx + 1)
y_coords = np.linspace(0, Ly, ny + 1)

nnx = nx + 1
nny = ny + 1
num_nodes = nnx * nny
num_elements = nx * ny

# Nodes
nodes = []
nid = 1
for j in range(nny):
    for i in range(nnx):
        nodes.append((nid, x_coords[i], y_coords[j]))
        nid += 1

# Elements
elements = []
eid = 1
for j in range(ny):
    for i in range(nx):
        n1 = j * nnx + i + 1
        n2 = n1 + 1
        n3 = n2 + nnx
        n4 = n1 + nnx
        elements.append((eid, n1, n2, n3, n4))
        eid += 1

# Node sets
# 1: left edge (x=0) — fix x-displacement
left_nodes = [j * nnx + 1 for j in range(nny)]
# 2: right edge (x=Lx) — applied displacement
right_nodes = [j * nnx + nnx for j in range(nny)]
# 3: bottom-left corner — fix y-displacement
bl_node = [1]
# 4: center line nodes (x=Lx/2, all y) — for phase-field d=1
center_i = nx // 2
center_nodes = [j * nnx + center_i + 1 for j in range(nny)]
# 5: all boundary nodes (for phase-field d=0 at far left/right)
far_pf_nodes = sorted(set(left_nodes + right_nodes))

node_sets = [
    (1, left_nodes),
    (2, right_nodes),
    (3, bl_node),
    (4, center_nodes),
    (5, far_pf_nodes),
]

with open("bar_coupled.geom", "w") as f:
    f.write("*version\n1.0\n")
    f.write("*title\n")
    f.write(f"Coupled test mesh: {nx}x{ny} quad4, Lx={Lx}, Ly={Ly}\n")
    f.write("*dimensions\n")
    f.write(f"{num_nodes}  # number of nodes\n")
    f.write(f"2   # number of spatial dimensions\n")
    f.write(f"1   # number of element sets\n")
    f.write(f"# [ID] [nel] [nen]\n")
    f.write(f"1   {num_elements}   4\n")
    f.write(f"{len(node_sets)}   # number of node sets\n")
    f.write(f"# [ID] [nnd]\n")
    for sid, nds in node_sets:
        f.write(f"{sid}   {len(nds)}\n")
    f.write(f"0   # number of side sets\n")
    f.write(f"# end dimensions\n")

    f.write("*nodesets\n")
    for sid, nds in node_sets:
        f.write("*set\n")
        f.write(f"{len(nds)}   # number of nodes\n")
        # Write in lines of 10
        for k in range(0, len(nds), 10):
            line = "  ".join(str(n) for n in nds[k:k+10])
            f.write(f"{line}\n")
    f.write("# end node sets\n")

    f.write("*sidesets\n")
    f.write("# end side sets\n")

    f.write("*elements\n")
    f.write("*set\n")
    f.write(f"{num_elements}  # number of elements\n")
    f.write(f"4  # number of element nodes\n")
    for e in elements:
        f.write(f"{e[0]}  {e[1]}  {e[2]}  {e[3]}  {e[4]}\n")
    f.write("# end elements\n")

    f.write("*nodes\n")
    f.write(f"{num_nodes}  # number of nodes\n")
    f.write(f"2   # number of spatial dimensions\n")
    for n in nodes:
        f.write(f"    {n[0]}   {n[1]:.10e}   {n[2]:.10e}\n")

print(f"Generated bar_coupled.geom:")
print(f"  Nodes: {num_nodes}, Elements: {num_elements}")
print(f"  NS 1 (left, {len(left_nodes)} nodes): fix D_X")
print(f"  NS 2 (right, {len(right_nodes)} nodes): applied D_X")
print(f"  NS 3 (bottom-left, 1 node): fix D_Y")
print(f"  NS 4 (center, {len(center_nodes)} nodes): d=1")
print(f"  NS 5 (far L/R, {len(far_pf_nodes)} nodes): d=0")
