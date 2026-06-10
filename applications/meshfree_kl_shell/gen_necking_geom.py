#!/usr/bin/env python3
"""Necking cylinder (paper Fig 12): R=10 (mid-surface), L=50, t=1. Stretched axially (z).
nt(circumferential, even) x nz(axial). Sets: 1=z=0 ring, 2=z=L ring, 9=all.
Usage: gen_necking_geom.py NT NZ"""
import sys, math
R, L = 10.0, 50.0
nt = int(sys.argv[1]) if len(sys.argv)>1 else 50
if nt%2: nt+=1
nz = int(sys.argv[2]) if len(sys.argv)>2 else 40
N = nt*nz
def nid(i,j): return j*nt+i+1
coords=[]
for j in range(nz):
    z=L*j/(nz-1)
    for i in range(nt):
        th=2*math.pi*i/nt; coords.append((nid(i,j), R*math.cos(th), R*math.sin(th), z))
elems=[]; e=1
for j in range(nz-1):
    for i in range(nt):
        i2=(i+1)%nt; elems.append((e,nid(i,j),nid(i2,j),nid(i2,j+1),nid(i,j+1))); e+=1
s1=[nid(i,0) for i in range(nt)]; s2=[nid(i,nz-1) for i in range(nt)]
s9=[nid(i,j) for j in range(nz) for i in range(nt)]
sets=[(1,s1),(2,s2),(9,s9)]
w=sys.stdout.write
w("*version\n1.0\n\n*title\nnecking %dx%d\n\n"%(nt,nz))
w("*dimensions\n%d\n3\n1\n1  %d  4\n%d\n"%(N,len(elems),len(sets)))
for sid,ns in sets: w("%d  %d\n"%(sid,len(ns)))
w("0\n\n*nodesets\n")
for sid,ns in sets: w("*set\n%d\n"%len(ns)); w("  ".join(map(str,ns))+"\n")
w("\n*sidesets\n\n*elements\n*set\n%d\n4\n"%len(elems))
for el in elems: w("  %d   %d %d %d %d\n"%el)
w("\n*nodes\n%d\n3\n"%N)
for (v,x,y,z) in coords: w("  %d   % .8e   % .8e   % .8e\n"%(v,x,y,z))
