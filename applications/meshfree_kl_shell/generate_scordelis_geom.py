#!/usr/bin/env python3
"""Generate a Tahoe .geom for the Scordelis-Lo roof as a SURFACE node cloud.

Cylindrical roof section: R=25, length L=50, 80 deg arc (+/-40 from crown).
Structured nt x nz grid of nodes on the surface (3D coords), with a background quad
mesh (for Tahoe's integration cells) and node sets for the diaphragm ends + the
free-edge midpoint (deflection probe). Membrane-dominated; reference free-edge
mid-span vertical deflection = 0.3006.

Node sets:
  1 = z=0 diaphragm ring
  2 = z=L diaphragm ring
  3 = single node for the axial (u_z) restraint (one corner)
  4 = free-edge midspan node (phi=+40, z=L/2) -- deflection probe
  9 = all nodes
Usage: python3 generate_scordelis_geom.py NT NZ > scordelis.geom
"""
import sys, math

nt = int(sys.argv[1]) if len(sys.argv) > 1 else 11
nz = int(sys.argv[2]) if len(sys.argv) > 2 else 11

R, L, phiMax = 25.0, 50.0, math.radians(40.0)
N = nt * nz
def nid(i, j):  # 1-based node id, i around arc, j along length
    return j * nt + i + 1

# coordinates (3D surface)
coords = []
for j in range(nz):
    for i in range(nt):
        phi = -phiMax + 2.0 * phiMax * i / (nt - 1)
        x = R * math.sin(phi)
        y = R * math.cos(phi)   # crown at y=R
        z = L * j / (nz - 1)
        coords.append((nid(i, j), x, y, z))

# background quad cells (for integration), 4-node, CCW
elems = []
eid = 1
for j in range(nz - 1):
    for i in range(nt - 1):
        elems.append((eid, nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)))
        eid += 1

# node sets
ns1 = [nid(i, 0)      for i in range(nt)]          # z=0 diaphragm
ns2 = [nid(i, nz - 1) for i in range(nt)]          # z=L diaphragm
ns3 = [nid(0, 0)]                                   # single axial restraint
ns4 = [nid(nt - 1, (nz - 1) // 2)]                 # free-edge midspan probe (phi=+40, z=L/2)
ns9 = [nid(i, j) for j in range(nz) for i in range(nt)]
nodesets = [(1, ns1), (2, ns2), (3, ns3), (4, ns4), (9, ns9)]

w = sys.stdout.write
w("*version\n1.0\n\n")
w("*title\nScordelis-Lo roof surface cloud %dx%d\n\n" % (nt, nz))
w("*dimensions\n")
w("%d\n3\n" % N)                                    # nodes, 3 spatial dims
w("1\n")                                            # 1 element set
w("1  %d  4\n" % len(elems))                        # set ID=1, nel, nen=4
w("%d\n" % len(nodesets))
for sid, ns in nodesets:
    w("%d  %d\n" % (sid, len(ns)))
w("0\n\n")                                          # 0 side sets

w("*nodesets\n")
for sid, ns in nodesets:
    w("*set\n%d\n" % len(ns))
    w("  ".join(str(n) for n in ns) + "\n")
w("\n*sidesets\n\n")

w("*elements\n*set\n%d\n4\n" % len(elems))
for e in elems:
    w("  %d   %d %d %d %d\n" % e)
w("\n*nodes\n%d\n3\n" % N)
for (nidv, x, y, z) in coords:
    w("  %d   % .8e   % .8e   % .8e\n" % (nidv, x, y, z))
