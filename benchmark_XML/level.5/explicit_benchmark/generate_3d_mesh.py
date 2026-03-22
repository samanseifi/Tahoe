#!/usr/bin/env python3
"""Generate 3D hexahedral meshes for explicit solver benchmarking.

Creates Tahoe .geom format meshes with varying element counts:
  - small:  10x5x5   = 250 elements
  - medium: 20x10x10 = 2000 elements
  - large:  40x20x10 = 8000 elements

Domain: [0, 2] x [0, 1] x [0, 1]
Single element block, 6 node sets (x=0, x=L, y=0, y=L, z=0, z=L).
"""
import os

def write_mesh(name, nx, ny, nz, Lx=2.0, Ly=1.0, Lz=1.0):
    nnx, nny, nnz = nx+1, ny+1, nz+1
    numnod = nnx * nny * nnz
    numel = nx * ny * nz

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)

    def node_id(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    # Node coordinates
    with open(f"{base}.geom.nd", "w") as f:
        f.write(f"{numnod}\n3\n")
        for k in range(nnz):
            for j in range(nny):
                for i in range(nnx):
                    nid = node_id(i, j, k)
                    x = i * Lx / nx
                    y = j * Ly / ny
                    z = k * Lz / nz
                    f.write(f"{nid:8d}   {x:.15e}   {y:.15e}   {z:.15e}\n")

    # Element connectivity (1-based, Hex8 standard ordering)
    with open(f"{base}.geom.es0", "w") as f:
        f.write(f"{numel}\n8\n")
        eid = 1
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    n1 = node_id(i,   j,   k)
                    n2 = node_id(i+1, j,   k)
                    n3 = node_id(i+1, j+1, k)
                    n4 = node_id(i,   j+1, k)
                    n5 = node_id(i,   j,   k+1)
                    n6 = node_id(i+1, j,   k+1)
                    n7 = node_id(i+1, j+1, k+1)
                    n8 = node_id(i,   j+1, k+1)
                    f.write(f"{eid:8d} {n1:7d} {n2:7d} {n3:7d} {n4:7d}"
                            f" {n5:7d} {n6:7d} {n7:7d} {n8:7d}\n")
                    eid += 1

    # Node sets: xmin(1), xmax(2), ymin(3), ymax(4), zmin(5), zmax(6)
    ns = [[], [], [], [], [], []]
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                nid = node_id(i, j, k)
                if i == 0:   ns[0].append(nid)  # x=0
                if i == nx:  ns[1].append(nid)  # x=Lx
                if j == 0:   ns[2].append(nid)  # y=0
                if j == ny:  ns[3].append(nid)  # y=Ly
                if k == 0:   ns[4].append(nid)  # z=0
                if k == nz:  ns[5].append(nid)  # z=Lz

    for idx in range(6):
        with open(f"{base}.geom.ns{idx}", "w") as f:
            f.write(f"{len(ns[idx])}\n")
            f.write("".join(f"{n:8d}" for n in ns[idx]) + "\n")

    # Side sets (one per face, referencing element-local face IDs)
    # Hex8 face numbering: 0=bottom(z-), 1=right(x+), 2=top(z+), 3=left(x-), 4=front(y-), 5=back(y+)
    # We use LS-DYNA/Tahoe convention for face ordering
    ss = [[], [], [], [], [], []]
    eid = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if i == 0:    ss[3].append((eid, 3))  # x=0 face
                if i == nx-1: ss[1].append((eid, 1))  # x=Lx face
                if j == 0:    ss[4].append((eid, 4))  # y=0 face
                if j == ny-1: ss[5].append((eid, 5))  # y=Ly face
                if k == 0:    ss[0].append((eid, 0))  # z=0 face
                if k == nz-1: ss[2].append((eid, 2))  # z=Lz face
                eid += 1

    for idx in range(6):
        with open(f"{base}.geom.ss{idx}", "w") as f:
            f.write(f"{len(ss[idx])}\n")
            for elem, face in ss[idx]:
                f.write(f"{elem:8d} {face:7d}\n")

    # Main .geom file
    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\nexplicit_3d_benchmark\n*dimensions\n")
        f.write(f"{numnod}\n3\n")
        f.write("1 # number of element sets\n")
        f.write(f"       1 {numel:8d}       8\n")
        f.write("6 # number of node sets\n")
        for idx in range(6):
            f.write(f"       {idx+1} {len(ns[idx]):7d}\n")
        f.write("6 # number of side sets\n")
        for idx in range(6):
            f.write(f"       {idx+1}       1 {len(ss[idx]):7d}\n")
        f.write("*nodesets\n")
        for idx in range(6):
            f.write(f"*set\n{name}.geom.ns{idx}\n")
        f.write("*sidesets\n")
        for idx in range(6):
            f.write(f"*set\n{name}.geom.ss{idx}\n")
        f.write("*elements\n*set\n")
        f.write(f"{name}.geom.es0\n")
        f.write("*nodes\n")
        f.write(f"{name}.geom.nd\n")

    print(f"{name}: {nx}x{ny}x{nz} = {numel} elements, {numnod} nodes")

if __name__ == "__main__":
    write_mesh("hex_small",  10, 5, 5)
    write_mesh("hex_medium", 20, 10, 10)
    write_mesh("hex_large",  40, 20, 10)
