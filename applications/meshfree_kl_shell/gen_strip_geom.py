#!/usr/bin/env python3
"""Flat strip in z=0 for a uniaxial tension test. Nx x Ny nodes, spacing 1.
Sets: 1=left edge(x=0), 2=right edge(x=Nx-1), 3=one left corner (rigid-y), 9=all."""
import sys
Nx=int(sys.argv[1]) if len(sys.argv)>1 else 21
Ny=int(sys.argv[2]) if len(sys.argv)>2 else 7
def nid(i,j): return j*Nx+i+1
coords=[(nid(i,j),float(i),float(j),0.0) for j in range(Ny) for i in range(Nx)]
elems=[]; e=1
for j in range(Ny-1):
    for i in range(Nx-1):
        elems.append((e,nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1))); e+=1
left =[nid(0,j) for j in range(Ny)]
right=[nid(Nx-1,j) for j in range(Ny)]
corner=[nid(0,0)]
alln=[nid(i,j) for j in range(Ny) for i in range(Nx)]
sets=[(1,left),(2,right),(3,corner),(9,alln)]
w=sys.stdout.write
w("*version\n1.0\n\n*title\nstrip %dx%d\n\n"%(Nx,Ny))
w("*dimensions\n%d\n3\n1\n1  %d  4\n%d\n"%(Nx*Ny,len(elems),len(sets)))
for sid,ns in sets: w("%d  %d\n"%(sid,len(ns)))
w("0\n\n*nodesets\n")
for sid,ns in sets: w("*set\n%d\n"%len(ns)); w("  ".join(map(str,ns))+"\n")
w("\n*sidesets\n\n*elements\n*set\n%d\n4\n"%len(elems))
for el in elems: w("  %d   %d %d %d %d\n"%el)
w("\n*nodes\n%d\n3\n"%(Nx*Ny))
for (v,x,y,z) in coords: w("  %d   % .8e   % .8e   % .8e\n"%(v,x,y,z))
