#!/usr/bin/env python3
"""Flat NxN node patch in the z=0 plane, spacing 1.0, with background quad cells.
For the stabilization unit test: a controlled patch where rigid/linear/hourglass modes are exact."""
import sys
N = int(sys.argv[1]) if len(sys.argv)>1 else 11
def nid(i,j): return j*N+i+1
coords=[(nid(i,j), float(i), float(j), 0.0) for j in range(N) for i in range(N)]
elems=[]; e=1
for j in range(N-1):
    for i in range(N-1):
        elems.append((e, nid(i,j), nid(i+1,j), nid(i+1,j+1), nid(i,j+1))); e+=1
allnodes=[nid(i,j) for j in range(N) for i in range(N)]
w=sys.stdout.write
w("*version\n1.0\n\n*title\nflat patch %dx%d\n\n"%(N,N))
w("*dimensions\n%d\n3\n1\n1  %d  4\n1\n9  %d\n0\n\n"%(N*N,len(elems),len(allnodes)))
w("*nodesets\n*set\n%d\n"%len(allnodes)); w("  ".join(str(n) for n in allnodes)+"\n")
w("\n*sidesets\n\n*elements\n*set\n%d\n4\n"%len(elems))
for el in elems: w("  %d   %d %d %d %d\n"%el)
w("\n*nodes\n%d\n3\n"%(N*N))
for (v,x,y,z) in coords: w("  %d   % .8e   % .8e   % .8e\n"%(v,x,y,z))
