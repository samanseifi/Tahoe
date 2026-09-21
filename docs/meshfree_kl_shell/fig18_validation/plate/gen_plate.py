#!/usr/bin/env python3
"""Square plate L x L in z=0, N x N nodes. Sets: 1=all edge nodes, 2=center node, 9=all."""
import sys
N=int(sys.argv[1]); L=float(sys.argv[2])
def nid(i,j): return j*N+i+1
h=L/(N-1)
coords=[(nid(i,j),i*h,j*h,0.0) for j in range(N) for i in range(N)]
elems=[]; e=1
for j in range(N-1):
    for i in range(N-1):
        elems.append((e,nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1))); e+=1
edge=[nid(i,j) for j in range(N) for i in range(N) if i in (0,N-1) or j in (0,N-1)]
center=[nid((N-1)//2,(N-1)//2)]
alln=[nid(i,j) for j in range(N) for i in range(N)]
sets=[(1,edge),(2,center),(9,alln)]
w=sys.stdout.write
w("*version\n1.0\n\n*title\nplate %dx%d L=%g\n\n"%(N,N,L))
w("*dimensions\n%d\n3\n1\n1  %d  4\n%d\n"%(N*N,len(elems),len(sets)))
for sid,ns in sets: w("%d  %d\n"%(sid,len(ns)))
w("0\n\n*nodesets\n")
for sid,ns in sets: w("*set\n%d\n"%len(ns)); w("  ".join(str(n) for n in ns)+"\n")
w("\n*sidesets\n\n*elements\n*set\n%d\n4\n"%len(elems))
for el in elems: w("  %d   %d %d %d %d\n"%el)
w("\n*nodes\n%d\n3\n"%(N*N))
for (v,x,y,z) in coords: w("  %d   % .8e   % .8e   % .8e\n"%(v,x,y,z))
