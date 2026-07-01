#!/usr/bin/env python3
"""Check that symmetric compression removed the loaded-edge bulge bias.
For the final frame of an exo:
  - top-surface rise should be mirror-symmetric about x=L/2  (left half ~ right half)
  - in-plane displacement DX should be ANTI-symmetric about x=L/2
We report the asymmetry of the top rise (max bulge left vs right) and the
left/right RMS of the rise; small + balanced => no biased bulge.
Usage: python3 check_symmetry.py <exo>
"""
import sys
import numpy as np
from netCDF4 import Dataset

f = sys.argv[1]
with Dataset(f) as ds:
    cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
    DX = np.array(ds["vals_nod_var1"][:]); DY = np.array(ds["vals_nod_var2"][:])
    t = np.array(ds["time_whole"][:])

top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
o = np.argsort(cx[top]); top = top[o]
xr = cx[top]; L = xr.max() - xr.min(); xc = xr.min() + L / 2.0
k = len(t) - 1                                  # final frame
y = cy[top] + DY[k, top]
rise = y - cy[top].mean()                        # vertical rise of the top surface
# mirror about center: compare rise(x) with rise(L-x)
xm = (xr - xc)
mirror = np.interp(-xm, xm, rise)                # rise at the mirrored location
asym = rise - mirror                             # 0 if perfectly symmetric

left = xm < 0; right = xm > 0
print(f"final frame t={t[k]:.1f}  L={L:.1f}  center x={xc:.1f}")
print(f"top rise:   mean={rise.mean():+.4f}  min={rise.min():+.4f}  max={rise.max():+.4f}")
print(f"left-half  rise RMS = {np.sqrt((rise[left]**2).mean()):.4f}")
print(f"right-half rise RMS = {np.sqrt((rise[right]**2).mean()):.4f}")
print(f"mirror asymmetry  RMS = {np.sqrt((asym**2).mean()):.4f}  (0 = perfectly symmetric)")
print(f"   -> asymmetry / rise-amplitude = "
      f"{np.sqrt((asym**2).mean())/max(np.sqrt((rise**2).mean()),1e-9):.2%}")
# peak bulge location (where is the material piling up?)
ipk = np.argmax(rise)
print(f"peak rise at x={xr[ipk]:.1f} ({(xr[ipk]-xr.min())/L:.0%} along film)"
      f"  [50% = center, biased if near 0% or 100%]")
