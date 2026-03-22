#!/usr/bin/env python3
"""Generate 2D quadrilateral meshes for explicit solver benchmarking.

Creates Tahoe .geom format meshes with varying element counts:
  - small:  20x10 = 200 elements
  - medium: 100x50 = 5000 elements
  - large:  200x100 = 20000 elements

Domain: [0, 2] x [0, 1]
Single element block, 4 node sets (bottom, right, top, left).
"""
import os

def write_mesh(name, nx, ny, Lx=2.0, Ly=1.0):
    """Write a single-block Q4 mesh in Tahoe .geom format."""
    nnx, nny = nx+1, ny+1
    numnod = nnx * nny
    numel = nx * ny

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)

    # Node coordinates
    with open(f"{base}.geom.nd", "w") as f:
        f.write(f"{numnod}\n2\n")
        for j in range(nny):
            for i in range(nnx):
                nid = j*nnx + i + 1
                x = i * Lx / nx
                y = j * Ly / ny
                f.write(f"{nid:8d}   {x:.15e}   {y:.15e}\n")

    # Element connectivity (1-based, counter-clockwise)
    with open(f"{base}.geom.es0", "w") as f:
        f.write(f"{numel}\n4\n")
        eid = 1
        for j in range(ny):
            for i in range(nx):
                n1 = j*nnx + i + 1
                n2 = n1 + 1
                n3 = n2 + nnx
                n4 = n1 + nnx
                f.write(f"{eid:8d} {n1:7d} {n2:7d} {n3:7d} {n4:7d}\n")
                eid += 1

    # Node sets: bottom (y=0), right (x=Lx), top (y=Ly), left (x=0)
    bottom = [j*nnx + i + 1 for i in range(nnx) for j in [0]]
    right  = [j*nnx + nx + 1 for j in range(nny)]
    top    = [ny*nnx + i + 1 for i in range(nnx)]
    left   = [j*nnx + 1 for j in range(nny)]

    for idx, nodes in enumerate([bottom, right, top, left]):
        with open(f"{base}.geom.ns{idx}", "w") as f:
            f.write(f"{len(nodes)}\n")
            f.write("".join(f"{n:8d}" for n in nodes) + "\n")

    # Side sets: bottom (face 0), right (face 1), top (face 2), left (face 3)
    # Format: element_local_id  face_id
    ss_bottom = [(i+1, 0) for i in range(nx)]  # bottom row, face 0
    ss_right  = [(j*nx + nx, 1) for j in range(ny)]  # right column, face 1
    ss_top    = [((ny-1)*nx + i + 1, 2) for i in range(nx)]  # top row, face 2
    ss_left   = [(j*nx + 1, 3) for j in range(ny)]  # left column, face 3

    for idx, sides in enumerate([ss_bottom, ss_right, ss_top, ss_left]):
        with open(f"{base}.geom.ss{idx}", "w") as f:
            f.write(f"{len(sides)}\n")
            for elem, face in sides:
                f.write(f"{elem:8d} {face:7d}\n")

    # Main .geom file
    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\nexplicit_benchmark\n*dimensions\n")
        f.write(f"{numnod}\n2\n")
        f.write("1 # number of element sets\n")
        f.write(f"       1 {numel:8d}       4\n")
        f.write("4 # number of node sets\n")
        for idx, nodes in enumerate([bottom, right, top, left]):
            f.write(f"       {idx+1} {len(nodes):7d}\n")
        f.write("4 # number of side sets\n")
        for idx, sides in enumerate([ss_bottom, ss_right, ss_top, ss_left]):
            f.write(f"       {idx+1}       1 {len(sides):7d}\n")
        f.write("*nodesets\n")
        for idx in range(4):
            f.write(f"*set\n{name}.geom.ns{idx}\n")
        f.write("*sidesets\n")
        for idx in range(4):
            f.write(f"*set\n{name}.geom.ss{idx}\n")
        f.write("*elements\n*set\n")
        f.write(f"{name}.geom.es0\n")
        f.write("*nodes\n")
        f.write(f"{name}.geom.nd\n")

    print(f"{name}: {nx}x{ny} = {numel} elements, {numnod} nodes")

if __name__ == "__main__":
    write_mesh("block_small",  20, 10)
    write_mesh("block_medium", 100, 50)
    write_mesh("block_large",  200, 100)
