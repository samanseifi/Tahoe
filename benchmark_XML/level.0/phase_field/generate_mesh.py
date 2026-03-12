#!/usr/bin/env python3
"""Generate a fine 1D-like quad4 mesh for the phase-field screened Poisson test.

Domain: x in [-L, L], y in [0, h] (thin strip)
The mesh is refined near x=0 where the solution has sharp gradients.

Boundary conditions:
  - d=1 at the center nodes (x=0)  -> Dirichlet
  - d=0 at the far ends (x=+-L)    -> Dirichlet (optional, natural BC suffices)

Analytical solution (1D): d(x) = exp(-|x|/ell)
"""

import numpy as np

# Parameters
L = 5.0       # half-length of domain
h = 0.1       # height (thin strip for 1D-like behavior)
nx = 100      # number of elements in x-direction
ny = 1        # single element in y-direction
ell = 1.0     # length scale (must match XML input)

# Generate x-coordinates: uniform spacing
x_coords = np.linspace(-L, L, nx + 1)

# y-coordinates
y_coords = np.linspace(0, h, ny + 1)

# Total nodes
nnx = nx + 1
nny = ny + 1
num_nodes = nnx * nny
num_elements = nx * ny

# Node coordinates (row-major: y varies slowest)
nodes = []
node_id = 1
for j in range(nny):
    for i in range(nnx):
        nodes.append((node_id, x_coords[i], y_coords[j]))
        node_id += 1

# Elements (4-node quads)
elements = []
elem_id = 1
for j in range(ny):
    for i in range(nx):
        n1 = j * nnx + i + 1
        n2 = n1 + 1
        n3 = n2 + nnx
        n4 = n1 + nnx
        elements.append((elem_id, n1, n2, n3, n4))
        elem_id += 1

# Node sets
# 1: bottom row (y=0)
bottom_nodes = list(range(1, nnx + 1))
# 2: top row (y=h)
top_nodes = list(range(nnx + 1, 2 * nnx + 1))
# 3: center nodes (x=0) -- for Dirichlet d=1
center_i = nx // 2  # index of x=0 node
center_nodes = [center_i + 1, center_i + 1 + nnx]
# 4: left edge (x=-L)
left_nodes = [1, nnx + 1]
# 5: right edge (x=+L)
right_nodes = [nnx, 2 * nnx]
# 6: all boundary nodes (left + right) for d=0
far_nodes = sorted(set(left_nodes + right_nodes))

node_sets = [
    (1, bottom_nodes),
    (2, top_nodes),
    (3, center_nodes),
    (4, left_nodes),
    (5, right_nodes),
    (6, far_nodes),
]

# Write geometry file
with open("bar_phase_field.geom", "w") as f:
    f.write("*version\n1.0\n")
    f.write("*title\n")
    f.write(f"Phase-field test mesh: {nx}x{ny} quad4, L={L}, h={h}\n")
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
        line = "  ".join(str(n) for n in nds)
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

print(f"Generated bar_phase_field.geom:")
print(f"  Nodes: {num_nodes}")
print(f"  Elements: {num_elements}")
print(f"  Center node set (ID=3): nodes {center_nodes}")
print(f"  Far boundary node set (ID=6): nodes {far_nodes}")
print(f"\nFor XML input:")
print(f"  Center BC: node_ID='3' (node set) -> d=1")
print(f"  Far BC:    node_ID='6' (node set) -> d=0")
print(f"  Or just use center BC with natural BC at far ends")
