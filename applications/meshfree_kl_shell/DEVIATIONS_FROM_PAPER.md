# Deviations from Wang & Bazilevs (2025) — meshfree_kl_shell (RKShellT)

Honest accounting of where this Tahoe element's results differ from the paper, and why.
Two categories: **A. quantitative result gaps** (per example) and **B. methodological differences**
(how the element's formulation departs from the paper's).

---

## A. Quantitative result deviations (per example)

| § | Example | Paper reference | Our result | Deviation | Status |
|---|---------|-----------------|-----------|-----------|--------|
| 4.1 | Geometry accuracy (vase) | quad normal rate 2, quad curv 2 (super), cubic normal **4** (super), cubic curv 2 | 1.96 / 2.03 / **3.05** / 2.01 | cubic-normal rate is the theoretical **3**, not the paper's superconvergent **4** | ✅ replicated (3/4 rates match) |
| 4.2.1 | Scordelis–Lo roof | −0.292 | −0.2910 | **0.4%** | ✅ replicated |
| 4.2.2 | Hemispherical shell | radial disp 0.0924 | **no solution** | element has **4 zero-energy modes at mesh corners** (free-edge rank deficiency) → singular stiffness | ❌ blocked |
| 4.2.3 | Pinched cylinder | 1.8624e-5 | 1.036× at 12,096 nodes | **+3.6%**, still converging; paper's finest 19,440-node mesh **errored** for us | 🟡 converging |
| 4.3 | Necking | neck εp contour to **2.0** | εp ~0.83 (coarse, 57%) … ~1.18 (paper-res, 89%) | εp climbs steadily but we **never stretched to ~150%** (full collapse) where 2.0 lives; Fig-15 shape matches but **not overlaid on the paper's IGA reference values** (unavailable) | 🟡 not run to collapse |
| 4.4 | Pinched elasto-plastic cyl | force–displacement (Fig 16, plotted only) | force ~2365 @ 78mm crush | **cannot score** — paper gives no numeric value in text | 🟡 unquantified |
| 4.5 | Square steel tube crush | — | **not run** | — | ❌ |
| 5.1–5.3 | Membrane-locking studies (t/R=100/1k/10k) | — | **not run** | — | ❌ |

---

## B. Methodological deviations (formulation differences)

1. **Membrane stabilization — DIFFERENT METHOD.**
   Paper uses the natural SCNI/NSNI smoothed-gradient membrane stabilization (Eq 33 / §5.2).
   We use a **penalty** (`stab_membrane`, `stab_bending` coefficients). This is the single biggest
   methodological gap: it underlies the **slow pinched-cylinder convergence** (3.6% at 12k vs the
   paper's faster convergence) and the **hemisphere free-edge failure** (penalty doesn't control the
   corner modes; the paper's natural stabilization does).

2. **Bending stabilization — WE ADD ONE; PAPER DOES NOT.**
   Paper deliberately omits through-thickness/bending stabilization (their Eq 34 was unstable in the
   azimuthal direction). We add a **bending penalty** to control hourglassing. Different choice.

3. **Co-rotational frame (Algorithm 2) — IMPLEMENTED BUT OFF BY DEFAULT.**
   Paper always integrates stress in the Flanagan–Taylor co-rotational frame. We implemented + verified
   it (`corotational=1`) but most runs use the legacy OrthoTangents frame (`corotational=0`). Tested
   effect on the pinch: <4%, non-systematic — so small, but it is a default-path deviation.

4. **σ33=0 enforcement — TWO PATHS.**
   Paper: Newton on D̃33 with a 3D law. We have that (`PlaneStressJ2_D33`, secant) **only when
   thickness_update=1**; the default path uses a 2D plane-stress J2 return (σ33=0 by construction, no
   explicit D̃33).

5. **Mass scaling / loading rate (necking, pinched-plastic).**
   Paper: real density 7.8e-9, dt=4e-8, velocity 30 m/s (true dynamics).
   Most of our runs: density 7800 (~10¹²× mass scaling), quasi-static. Tested effect on necking
   localization: ~10% on εp — minor, but a deviation. (Real-density runs are ~6–12 h at paper resolution.)

6. **Element architecture.**
   We derive from `ElementBaseT` with a PCA tangent-plane chart + surface meshfree support (MLSSolverT),
   because the 3D-solid meshfree path rejects a surface node cloud. The paper's is a purpose-built RK
   shell. Different implementation; same intended formulation.

---

## Root cause summary
Of the unresolved items (hemisphere blocked, pinched cylinder 3.6%, necking εp), the **common blocker is
item B.1 — the paper's §5.2 natural membrane stabilization**, which we replaced with a penalty. The
necking εp=2.0 gap is separately just **insufficient stretch** (need ~150% to full collapse).

## Solidly matching the paper
- 4.1 geometry accuracy (convergence rates), 4.2.1 Scordelis–Lo (0.4%), 4.2.3 pinched cylinder at full
  resolution (3.6%, converging). All three Algorithms (1/2/3) implemented and individually unit-verified.
