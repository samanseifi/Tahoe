#!/usr/bin/env python3
"""Vase shape (paper sec 4.1, Eq 100): rho(z)=R+5 sin(2 pi z/L), R=25, L=50.
x=(rho cos th, rho sin th, z). nt(circumferential, even) x nz(axial). Set 9=all.
Usage: gen_vase_geom.py NT NZ"""
import sys, math
R, L = 25.0, 50.0
nt = int(sys.argv[1]) if len(sys.argv)>1 else 40
if nt%2: nt+=1
nz = int(sys.argv[2]) if len(sys.argv)>2 else 40
N = nt*nz
def nid(i,j): return j*nt+i+1
def rho(z): return R + 5.0*math.sin(2.0*math.pi*z/L)
coords=[]
for j in range(nz):
    z=L*j/(nz-1); rr=rho(z)
    for i in range(nt):
        th=2*math.pi*i/nt; coords.append((nid(i,j), rr*math.cos(th), rr*math.sin(th), z))
elems=[]; e=1
for j in range(nz-1):
    for i in range(nt):
        i2=(i+1)%nt; elems.append((e,nid(i,j),nid(i2,j),nid(i2,j+1),nid(i,j+1))); e+=1
s9=[nid(i,j) for j in range(nz) for i in range(nt)]
sets=[(9,s9)]
w=sys.stdout.write
w("*version\n1.0\n\n*title\nvase %dx%d\n\n"%(nt,nz))
w("*dimensions\n%d\n3\n1\n1  %d  4\n%d\n"%(N,len(elems),len(sets)))
for sid,ns in sets: w("%d  %d\n"%(sid,len(ns)))
w("0\n\n*nodesets\n")
for sid,ns in sets: w("*set\n%d\n"%len(ns)); w("  ".join(map(str,ns))+"\n")
w("\n*sidesets\n\n*elements\n*set\n%d\n4\n"%len(elems))
for el in elems: w("  %d   %d %d %d %d\n"%el)
w("\n*nodes\n%d\n3\n"%N)
for (v,x,y,z) in coords: w("  %d   % .8e   % .8e   % .8e\n"%(v,x,y,z))
