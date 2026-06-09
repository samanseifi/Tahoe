# Fig 18 (pinched elasto-plastic cylinder, Wang & Bazilevs §4.4) — status & open problem

## What is reproduced (solid)
- **Material calibrated to the paper**: E=3.0e3 MPa, ν=0.3, linear isotropic hardening
  K(εp)=24.3+300·εp. (My deck had all three values 1000× too high — fixed.)
- **Geometry/BC**: 8320 nodes (nt=160×nz=52), L=600, R=300, h=3, fixed diaphragm ends,
  two diametric mid-length **pinch nodes** (the paper's load — NOT a line load).
- **Qualitative behavior correct**: cross-section folds into the butterfly/bowtie; the
  EQ_PLASTIC_STRAIN field (clip 0–0.1) shows the distributed plastic-hinge pattern of Fig 17
  (99% yielded, max εp≈0.43); load-displacement has the rise→peak→post-buckling-drop→densification
  shape.

## Fixes applied along the way (all committed)
1. Bending-hourglass control (curvature residual) → arrested the inextensional sawtooth.
2. Boundary z-anchor moved off the diaphragm ring to the load points (a single-node z constraint on
   a meshfree boundary ring hourglassed that ring; u_x=u_y=0 verified exact at both ends).
3. Frequency-selective bending operator R = D_chord − D_LS (moment-matched LS Laplacian reference).
4. EQ_PLASTIC_STRAIN output (conditional on fYield>0; benchmark stays 3-field bit-exact).
5. Total-reaction sum over a load set (monitor_stride/monitor_count) for distributed loads.
6. include_bend flag → report the MATERIAL reaction (exclude the bending-penalty from the reported force).

## The open problem: force magnitude is ~4× too high
Paper Fig 18 ≈ <1000 at 150mm; this element ≈ 4000 (peak 4670 @ 204mm). **Diagnosed, not mysterious:**

- It is NOT a measurement bug. The reaction quantity is correct.
- It is the **bending-control over-stiffening the deformation.** Evidence (β-sweep @ ~36mm):
  | β   | material force | mid-ring smoothness |
  |-----|----------------|---------------------|
  | 0   | ~480 (paper scale) | 2.6 (hourglass)  |
  | 1   | 1356           | 1.11 (hourglass)    |
  | 2   | 1513           | 0.71                |
  | 8   | 1641 → peak 4670 | 0.10 (smooth)     |
  There is **no β sweet spot**: the moment β tames the hourglass, the force has already jumped 3×,
  and it's still hourglassing below that.

## Root cause
The bending hourglass is the standard nodal-integration deficiency (the RKPM kernel smooths the
node-to-node curvature → under-integrated → zero-energy mode; the paper has it too). My stabilization
is a **penalty ∝ β·R²**, which blows up on the **sub-grid sharp features** (the point-load tip and
the plastic hinges) — a curvature penalty cannot distinguish a one-node-wide physical feature from a
hourglass. The paper avoids this entirely: it uses **no bending curvature penalty**. Its 3
through-thickness points carry the bending stiffness, and its **consistent natural (Taylor)
stabilization** — scaling as (cell size)²/12 × strain-gradient, vanishing on smooth fields — removes
the hourglass *without* adding stiffness to the folding.

## The fix (a genuine reformulation, NOT a parameter tweak)
Replace the R² curvature penalty with a **consistent stabilization**. Two routes:
- **Full natural/Taylor stabilization (paper Eq. 33)**: the in-plane gradient of the strain energy
  × M_K=(cell size)²/12. For the *bending* part this needs **3rd-order RKPM derivatives** (the paper
  flags these as costly and "may lead to instability"); requires extending MLSSolverT to 3rd order.
- **Normal-gradient curvature (cheaper, no 3rd derivatives)**: build the curvature from the
  divergence-theorem of the **deformed nodal normals** (1st-derivative quantities that *do* oscillate
  on the sawtooth). Concrete plan:
  1. Two-pass RHS: pass 1 computes a deformed unit-normal field per node (n_K = x,1×x,2 from
     fXref+fDphi+u); pass 2 uses it.
  2. Normal-based curvature κ_n = div of the nodal-normal tilt; residual R = κ_n − κ_resolved.
  3. **Scale by the natural M_K=(cell size)²/12** (small, bounded) — NOT β·R² — so sharp sub-grid
     features stay bounded.
  4. Caveat: κ_n is **nonlinear in u**, so the force is coeff·R·dR/du (not the current rank-1 linear
     penalty); and it is uncertain whether even this fully resolves the sub-grid tension.

## Recommendation
Tackle the natural-stabilization reformulation as a **dedicated, focused effort** (the normal-gradient
route is the principled, cheaper one). Until then: the element reproduces §4.4 **qualitatively**
(correct material, butterfly mode, distributed hinges); the **force magnitude is over-stiffened ~4×
by the curvature-penalty stabilization — a documented, understood limitation, not a bug.**

## RESOLUTION (the over-stiffening is FIXED)
The 4x force over-stiffening is resolved by replacing the curvature PENALTY with the paper's actual
**Eq. 33 natural (Taylor-gradient) stabilization** -- which already existed as kernels
(BMatrixGradient + BMatrixCurvatureGradient in KLShellKernels.h, using the 3rd derivatives the MLS
already computes via DDDphi) but was never wired in. Now wired (membrane + bending Taylor, stab_natural
param). Because it is CONSISTENT (vanishes on smooth fields), it does NOT over-stiffen:
- GATE 1 (consistency): Scordelis-Lo bit-exact -0.292, KLShell 22/22, with the stabilizer ON.
- GATE 2 (force): per-load-point reaction ~432 (paper scale), NOT the 4670 penalty artifact.

One physical subtlety: a SINGLE-NODE pinch is a sub-grid delta that excites a kernel-ZERO-ENERGY
inextensional sawtooth -- the RKPM kernel smooths it to zero, so NO kernel-based stabilizer (the
paper's included) can catch it. A RESOLVED load (3x3 patch -> smoothness 0.34; or the theta=0/pi
generator lines = the paper's "top and bottom SURFACE" -> smoothness 0.16) does not excite it.
What IS resolved: the over-stiffening MECHANISM (the penalty inflation) is gone; the stabilizer is
consistent (Gate 1) and the early-crush force is paper-scale (432/node @24mm).

## Honest open item (full line curve, 7 pts to 156mm)
The LINE load is SMOOTH (0.017 @156mm) but it is the WRONG problem: it crushes the whole 600mm length
uniformly -> global ovalization -> per-node force rises MONOTONICALLY 432->2250 (24->156mm), crossing
the paper's <1000 at ~70mm and hitting ~2x at 150mm, with NO post-buckling drop. The paper's curve is
a LOCALIZED buckling (rise->peak->drop->~<1000). So the line proves smoothness+consistency but does
NOT reproduce Fig 18's shape or 150mm magnitude. The missing ingredient is the LOAD FOOTPRINT: a load
localized enough to BUCKLE (the drop) yet resolved enough not to sawtooth. Single node buckles but
sawtooths; line is smooth but won't buckle. The 3x3 PATCH is the in-between candidate -- under test
(does its force peak+drop while staying smooth?).
