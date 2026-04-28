#!/usr/bin/env python3
"""Generate Tet4 meshes by splitting each Hex8 cell into 5 tets.

Used to verify the explicit Tet4 kernel runs end-to-end and to seed
benchmarks for ANP/F-bar Tet4 (#28) — a hex8 mesh of size NxNxN gives
5*N^3 tets, each cell decomposed via the standard 5-tet pattern.

Decomposition (one of several possible):
  a hex with corner indices (000)..(111), labelled
       0=(0,0,0), 1=(1,0,0), 2=(1,1,0), 3=(0,1,0)   bottom
       4=(0,0,1), 5=(1,0,1), 6=(1,1,1), 7=(0,1,1)   top
  splits into 5 tets:
       (0, 1, 2, 5)
       (0, 2, 3, 7)
       (0, 5, 7, 4)
       (2, 5, 6, 7)
       (0, 2, 5, 7)        # central tet
"""
import os

# Standard 5-tet split of a hex8.  Each row is one tet.
# Hex8 corner ordering: bottom 0,1,2,3 then top 4,5,6,7 (Tahoe convention).
# Tet node order chosen so that J = (x_0-x_2, x_1-x_2, x_3-x_2) has positive
# determinant under Tahoe's Tet4 parametric layout (node 0 at r=1, node 1 at
# s=1, node 2 at origin, node 3 at t=1).  Verified algebraically:
#   tets 0..3 each have detJ=+1, central tet has detJ=+2; total = 6, matches
#   the unit hex volume of 1.
HEX_TO_TETS = [
    (0, 1, 2, 5),
    (0, 2, 3, 7),
    (0, 5, 7, 4),
    (2, 5, 6, 7),
    (0, 2, 7, 5),  # central — note 7,5 ordering for positive Jacobian
]


def write_tet_mesh(name, nx, ny, nz, Lx=1.0, Ly=1.0, Lz=1.0):
    nnx, nny, nnz = nx + 1, ny + 1, nz + 1
    numnod = nnx * nny * nnz
    numel_hex = nx * ny * nz
    numel_tet = numel_hex * len(HEX_TO_TETS)

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)

    def nid(i, j, k):
        return k * nnx * nny + j * nnx + i + 1

    # node coordinates
    with open(f"{base}.geom.nd", "w") as f:
        f.write(f"{numnod}\n3\n")
        for k in range(nnz):
            for j in range(nny):
                for i in range(nnx):
                    x = i * Lx / nx
                    y = j * Ly / ny
                    z = k * Lz / nz
                    f.write(f"{nid(i,j,k):8d}   {x:.15e}   {y:.15e}   {z:.15e}\n")

    # element connectivity (Tet4)
    with open(f"{base}.geom.es0", "w") as f:
        f.write(f"{numel_tet}\n4\n")
        eid = 1
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    # hex corner node IDs (Tahoe convention: bottom 0-3 ccw, top 4-7 ccw)
                    hex_n = [
                        nid(i,   j,   k),    # 0
                        nid(i+1, j,   k),    # 1
                        nid(i+1, j+1, k),    # 2
                        nid(i,   j+1, k),    # 3
                        nid(i,   j,   k+1),  # 4
                        nid(i+1, j,   k+1),  # 5
                        nid(i+1, j+1, k+1),  # 6
                        nid(i,   j+1, k+1),  # 7
                    ]
                    for tet in HEX_TO_TETS:
                        # Tahoe Tet4 node convention (per Tet4KernelT header):
                        # node 0 at r=1, node 1 at s=1, node 2 at origin, node 3 at t=1.
                        # Use the natural HEX_TO_TETS ordering — the Jacobian sign
                        # is handled by detJ (tet mesh decomp may produce mixed signs;
                        # Tet4KernelT uses fabs in volume but signed det in dN/dx).
                        n = [hex_n[v] for v in tet]
                        f.write(f"{eid:8d} {n[0]:7d} {n[1]:7d} {n[2]:7d} {n[3]:7d}\n")
                        eid += 1

    # node sets: 6 face sets + 1 body (everything not on z=0)
    ns = [[], [], [], [], [], [], []]
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                n = nid(i, j, k)
                if i == 0:   ns[0].append(n)
                if i == nx:  ns[1].append(n)
                if j == 0:   ns[2].append(n)
                if j == ny:  ns[3].append(n)
                if k == 0:   ns[4].append(n)
                if k == nz:  ns[5].append(n)
                if k > 0:    ns[6].append(n)   # body (non-impact-face)

    for idx in range(7):
        with open(f"{base}.geom.ns{idx}", "w") as f:
            f.write(f"{len(ns[idx])}\n")
            f.write("".join(f"{n:8d}" for n in ns[idx]) + "\n")

    # main .geom file
    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\ntet_explicit_benchmark\n*dimensions\n")
        f.write(f"{numnod}\n3\n")
        f.write("1 # number of element sets\n")
        f.write(f"       1 {numel_tet:8d}       4\n")
        f.write("7 # number of node sets\n")
        for idx in range(7):
            f.write(f"       {idx+1} {len(ns[idx]):7d}\n")
        f.write("0 # number of side sets\n")
        f.write("*nodesets\n")
        for idx in range(7):
            f.write(f"*set\n{name}.geom.ns{idx}\n")
        f.write("*sidesets\n")
        f.write("*elements\n*set\n")
        f.write(f"{name}.geom.es0\n")
        f.write("*nodes\n")
        f.write(f"{name}.geom.nd\n")

    print(f"{name}: {nx}x{ny}x{nz} hex => {numel_tet} tets, {numnod} nodes")


if __name__ == "__main__":
    write_tet_mesh("tet_small",  3, 3, 6)    # 3*3*6=54 hex => 270 tets
    write_tet_mesh("tet_medium", 5, 5, 10)   # 250 hex => 1250 tets
