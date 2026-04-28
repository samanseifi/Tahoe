#!/usr/bin/env python3
"""Generate Taylor bar hexahedral meshes (quarter-symmetry).

Classic Taylor bar impact benchmark, Wilkins & Guinan 1973:
  - Cu cylinder of length L0 = 32.4 mm, diameter D0 = 6.4 mm
  - Impact velocity V0 = 190 m/s against rigid wall
  - Reference final length Lf ~= 21.47 mm (33.8% compression)
  - Reference mushroom diameter ~= 13.5 mm (~2.1x initial)

We use a square-cross-section prism approximation with quarter-symmetry:
  - Domain: (x,y,z) in [0, a/2] x [0, a/2] x [0, L0]
  - a = sqrt(pi) * D0/2 preserves cross-sectional area  (equiv-square)
  - Symmetry BC: x=0 face has u_x=0, y=0 face has u_y=0
  - Impact face: z=0 has u_z=0 (rigid frictionless wall)

Node set IDs:
  1: x=0 face (symmetry: u_x=0)
  2: y=0 face (symmetry: u_y=0)
  3: z=0 face (rigid wall: u_z=0)
  4: all nodes (for initial velocity application)

Units convention: mm, ms, g, MPa (consistent: F=ma -> MPa * mm^2 = N, rho g/mm^3)
  - Length: mm, time: ms, mass: g
  - Velocity: mm/ms = m/s
  - Stress: g/(mm*ms^2) = MPa
  - Density: g/mm^3 -> Cu: 8.96e-3
"""
import os


def write_mesh(name: str, nr: int, nz: int,
               L0: float = 32.4, a: float = 5.67):
    """Generate quarter-symmetry prismatic Taylor bar mesh.

    nr: elements along each radial side (x and y same)
    nz: elements along axial (z) direction
    a:  full cross-section side length (mm).  Default: sqrt(pi)*D0/2 ~ 5.67 mm
        preserves area of D0=6.4mm round cross section.
    L0: bar length (mm).  Default 32.4 (Wilkins & Guinan 1973).
    """
    half_a = a / 2.0
    nnr, nnz = nr + 1, nz + 1
    numnod = nnr * nnr * nnz
    numel = nr * nr * nz

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)

    def nid(i, j, k):
        return k * nnr * nnr + j * nnr + i + 1

    # node coordinates
    with open(f"{base}.geom.nd", "w") as f:
        f.write(f"{numnod}\n3\n")
        for k in range(nnz):
            for j in range(nnr):
                for i in range(nnr):
                    x = i * half_a / nr
                    y = j * half_a / nr
                    z = k * L0 / nz
                    f.write(f"{nid(i,j,k):8d}   {x:.15e}   {y:.15e}   {z:.15e}\n")

    # element connectivity (Hex8)
    with open(f"{base}.geom.es0", "w") as f:
        f.write(f"{numel}\n8\n")
        eid = 1
        for k in range(nz):
            for j in range(nr):
                for i in range(nr):
                    n = [nid(i,   j,   k),   nid(i+1, j,   k),
                         nid(i+1, j+1, k),   nid(i,   j+1, k),
                         nid(i,   j,   k+1), nid(i+1, j,   k+1),
                         nid(i+1, j+1, k+1), nid(i,   j+1, k+1)]
                    f.write(f"{eid:8d} " + " ".join(f"{x:7d}" for x in n) + "\n")
                    eid += 1

    # node sets: 0=x=0, 1=y=0, 2=z=0, 3=z=L0, 4=all nodes NOT on z=0 (for initial velocity)
    ns = [[], [], [], [], []]
    for k in range(nnz):
        for j in range(nnr):
            for i in range(nnr):
                n = nid(i, j, k)
                if i == 0:  ns[0].append(n)  # x=0 symmetry
                if j == 0:  ns[1].append(n)  # y=0 symmetry
                if k == 0:  ns[2].append(n)  # z=0 impact face
                if k == nz: ns[3].append(n)  # z=L0 free end
                if k > 0:   ns[4].append(n)  # body (everything except wall)

    for idx in range(5):
        with open(f"{base}.geom.ns{idx}", "w") as f:
            f.write(f"{len(ns[idx])}\n")
            f.write("".join(f"{n:8d}" for n in ns[idx]) + "\n")

    # Master .geom file (inline references to sub-files)
    with open(f"{base}.geom", "w") as f:
        f.write("*version\n1.0\n*title\ntaylor_bar_quarter_symmetry\n*dimensions\n")
        f.write(f"{numnod}\n3\n")
        f.write("1 # number of element sets\n")
        f.write(f"       1 {numel:8d}       8\n")
        f.write("5 # number of node sets\n")
        for idx in range(5):
            f.write(f"       {idx+1} {len(ns[idx]):7d}\n")
        f.write("0 # number of side sets\n")
        f.write("*nodesets\n")
        for idx in range(5):
            f.write(f"*set\n{name}.geom.ns{idx}\n")
        f.write("*sidesets\n")
        f.write("*elements\n*set\n")
        f.write(f"{name}.geom.es0\n")
        f.write("*nodes\n")
        f.write(f"{name}.geom.nd\n")

    print(f"{name}: {nr}x{nr}x{nz} = {numel} Hex8 elems, {numnod} nodes"
          f" (cross {a:.2f}x{a:.2f} mm, length {L0:.1f} mm)")


if __name__ == "__main__":
    # Small mesh for fast smoke test / CI (runs in ~1-2 seconds)
    write_mesh("taylor_small",  nr=3,  nz=15)   # 3x3x15 =  135 elements
    # Medium mesh for benchmarking (~10-30 seconds)
    write_mesh("taylor_medium", nr=5,  nz=25)   # 5x5x25 =  625 elements
    # Larger mesh for validation (~1 minute)
    write_mesh("taylor_large",  nr=8,  nz=40)   # 8x8x40 = 2560 elements
