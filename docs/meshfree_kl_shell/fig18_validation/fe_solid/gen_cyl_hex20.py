#!/usr/bin/env python3
"""1/8-symmetry solid (hex20) mesh of the pinched cylinder: R=300 (mid-surface), h=3, L=600.
Octant: theta in [0,pi/2], z in [0,L/2], radial through the thickness. Pinch point at theta=0, z=L/2.
Node sets: 1 = z=0 diaphragm (u_x=u_y=0)   2 = z=L/2 symmetry (u_z=0)   3 = theta=0 plane (u_y=0)
           4 = theta=pi/2 plane (u_x=0)      5 = pinch line (theta=0, z=L/2, all through-thickness nodes)
Usage: gen_cyl_hex20.py NTHETA NZ NR > cyl.geom   (full-load = 4 x reaction of set 5)
Hex20 local ordering (Tahoe/Exodus): vertices bottom(t-) (r-,s-),(r+,s-),(r+,s+),(r-,s+), top(t+) same;
mid-edges 9-12 bottom ring, 13-16 top ring, 17-20 vertical. Here r<->theta, s<->z, t<->radial."""
import sys, math
R, h, L = 300.0, 3.0, 600.0
nt, nz, nr = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
pgrade = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0   # >1: grade spacing finer toward the pinch (theta=0, z=L/2)
Nt, Nz, Nr = 2*nt+1, 2*nz+1, 2*nr+1      # quadratic lattice
ids = {}; coords = []
def lat(i, j, k):
    odd = (i % 2) + (j % 2) + (k % 2)
    return odd <= 1                          # serendipity: skip face-center and body-center nodes
for k in range(Nr):
    for j in range(Nz):
        for i in range(Nt):
            if not lat(i, j, k): continue
            th = (math.pi/2) * (i/(Nt-1))**pgrade; z = (L/2) * (1.0 - (1.0 - j/(Nz-1))**pgrade); r = R - h/2 + h * k / (Nr-1)
            ids[(i, j, k)] = len(coords) + 1
            coords.append((r*math.cos(th), r*math.sin(th), z))
elems = []
for et in range(nt):
    for ez in range(nz):
        for er in range(nr):
            i0, j0, k0 = 2*et, 2*ez, 2*er
            v = lambda di, dj, dk: ids[(i0+di, j0+dj, k0+dk)]
            n = [v(0,0,0), v(2,0,0), v(2,2,0), v(0,2,0), v(0,0,2), v(2,0,2), v(2,2,2), v(0,2,2),
                 v(1,0,0), v(2,1,0), v(1,2,0), v(0,1,0),
                 v(1,0,2), v(2,1,2), v(1,2,2), v(0,1,2),
                 v(0,0,1), v(2,0,1), v(2,2,1), v(0,2,1)]
            elems.append(n)
def pick(cond): return [ids[key] for key in sorted(ids, key=lambda q: ids[q]) if cond(*key)]
s1 = pick(lambda i,j,k: j == 0)
s2 = pick(lambda i,j,k: j == Nz-1)
s3 = pick(lambda i,j,k: i == 0)
s4 = pick(lambda i,j,k: i == Nt-1)
s5 = pick(lambda i,j,k: i == 0 and j == Nz-1)
sets = [(1,s1),(2,s2),(3,s3),(4,s4),(5,s5)]
w = sys.stdout.write
w("*version\n1.0\n\n*title\npinched cylinder octant hex20 %dx%dx%d\n\n" % (nt, nz, nr))
w("*dimensions\n%d\n3\n1\n1  %d  20\n%d\n" % (len(coords), len(elems), len(sets)))
for sid, ns in sets: w("%d  %d\n" % (sid, len(ns)))
w("0\n\n*nodesets\n")
for sid, ns in sets:
    w("*set\n%d\n" % len(ns)); w("  ".join(str(n) for n in ns) + "\n")
w("\n*sidesets\n\n*elements\n*set\n%d\n20\n" % len(elems))
for e, n in enumerate(elems, 1): w("  %d   %s\n" % (e, " ".join(str(x) for x in n)))
w("\n*nodes\n%d\n3\n" % len(coords))
for v, (x, y, z) in enumerate(coords, 1): w("  %d   % .10e   % .10e   % .10e\n" % (v, x, y, z))
sys.stderr.write("nodes %d elements %d pinch-line nodes %d\n" % (len(coords), len(elems), len(s5)))
