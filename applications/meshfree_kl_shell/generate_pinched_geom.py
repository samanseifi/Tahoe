#!/usr/bin/env python3
"""Generate a Tahoe .geom for the pinched cylinder with end diaphragms (full model).

Belytschko obstacle course: R=300, L=600, t=3, E=3e6, nu=0.3, two diametric unit point
loads at mid-length. Reference radial deflection under the load = 1.8248e-5.

Structured nt (circumferential, EVEN) x nz (axial, ODD) node cloud on the cylinder surface,
background quad cells, and node sets:
  1 = z=0 diaphragm ring         2 = z=L diaphragm ring
  3 = load point A (theta=0,   z=L/2)   -> load in -x (inward)
  4 = load point B (theta=pi,  z=L/2)   -> load in +x (inward)
  5 = single node (theta=0, z=0) to remove the axial (z) rigid-body mode
  9 = all nodes
Usage: python3 generate_pinched_geom.py NT NZ > pinched.geom   (NT even, NZ odd)
"""
import sys, math

R, L = 300.0, 600.0

# nt = circumferential count (EVEN, node at theta=pi). nz auto-balances to ~square cells
# (arc spacing 2piR/nt ~ axial spacing L/(nz-1)) unless given explicitly. A balanced cloud is
# required: the isotropic kernel support must reach neighbors in BOTH directions, else an
# unbalanced mesh disconnects into axial strips (singular system).
import math as _m
nt = int(sys.argv[1]) if len(sys.argv) > 1 else 24
if nt % 2 != 0: nt += 1
if len(sys.argv) > 2:
    nz = int(sys.argv[2])
else:
    nz = int(round(nt * L / (2.0 * _m.pi * R))) + 1   # axial count for square cells
if nz % 2 == 0: nz += 1          # need a node at z=L/2
if nz < 5: nz = 5
N = nt * nz
def nid(i, j): return j * nt + i + 1     # 1-based; i around, j along

coords = []
for j in range(nz):
    z = L * j / (nz - 1)
    for i in range(nt):
        th = 2.0 * math.pi * i / nt
        coords.append((nid(i, j), R*math.cos(th), R*math.sin(th), z))

# background quad cells (wrap around in theta)
elems = []; eid = 1
for j in range(nz - 1):
    for i in range(nt):
        i2 = (i + 1) % nt
        elems.append((eid, nid(i, j), nid(i2, j), nid(i2, j+1), nid(i, j+1))); eid += 1

jmid = (nz - 1) // 2
ns1 = [nid(i, 0)      for i in range(nt)]
ns2 = [nid(i, nz - 1) for i in range(nt)]
ns3 = [nid(0,      jmid)]      # theta=0
ns4 = [nid(nt // 2, jmid)]     # theta=pi
ns5 = [nid(0, 0)]
ns9 = [nid(i, j) for j in range(nz) for i in range(nt)]
nodesets = [(1,ns1),(2,ns2),(3,ns3),(4,ns4),(5,ns5),(9,ns9)]

w = sys.stdout.write
w("*version\n1.0\n\n*title\nPinched cylinder %dx%d\n\n" % (nt, nz))
w("*dimensions\n%d\n3\n1\n1  %d  4\n%d\n" % (N, len(elems), len(nodesets)))
for sid, ns in nodesets: w("%d  %d\n" % (sid, len(ns)))
w("0\n\n*nodesets\n")
for sid, ns in nodesets:
    w("*set\n%d\n" % len(ns)); w("  ".join(str(n) for n in ns) + "\n")
w("\n*sidesets\n\n*elements\n*set\n%d\n4\n" % len(elems))
for e in elems: w("  %d   %d %d %d %d\n" % e)
w("\n*nodes\n%d\n3\n" % N)
for (v,x,y,z) in coords: w("  %d   % .8e   % .8e   % .8e\n" % (v,x,y,z))
