# Uniaxial tension — J2 cubic-spline hardening verification (#49)

Pulls a slim quarter-symmetry bar in tension to verify that the
`<Simo_J2>` radial-return reproduces the input `cubic_spline` σ_y(ε_p)
hardening curve under monotonic uniaxial loading.  Companion to the
Brinell benchmark (`../brinell/`) which uses the same material under
contact loading — this benchmark isolates the material from contact
nonlinearities.

## Files

| File | Notes |
|------|-------|
| [`generate_tensile_mesh.py`](generate_tensile_mesh.py) | Writes `tensile.geom`, a 4 × 4 × 40 Hex8 quarter bar (640 elements). |
| [`tensile.xml`](tensile.xml) | 25 quasi-static steps to δ = 0.5 mm (= 10 % engineering strain). |
| [`compare_to_hardening.py`](compare_to_hardening.py) | Extracts (α, σ_zz, σ_VM) history at the bar-centre probe and overlays the input spline. |

## Setup

- **Geometry** (quarter-symmetry): 0.5 × 0.5 × 5 mm; full bar would be 1 × 1 × 10 mm.
- **Mesh**: 4 × 4 × 40 Hex8 = 640 elements, 1025 nodes.
- **Material**: same family as `level.5/brinell/`, with a denser knot set
  in the rapidly-changing region of the spline:
  ```
  Simo_J2,  E = 200 000 MPa,  ν = 0.3,  ρ = 7.85
  cubic_spline σ_y(ε_p):
      (0.000, 250)  (0.005, 320)  (0.020, 420)  (0.050, 500)
      (0.080, 535)  (0.110, 560)  (0.150, 580)  (0.200, 600)
      (0.300, 625)  (0.400, 640)  (0.500, 650)         [MPa]
  ```
- **BCs**: top face (NS1) prescribed `u_z = 0.5 mm`, bottom face (NS2)
  `u_z = 0`, symmetry on x = 0 (NS3) and y = 0 (NS4).  Lateral faces
  free → Poisson contraction.
- **Solver**: `<nonlinear_solver_LS>` + SPOOLES, line search active.

## Running

```bash
cd benchmark_XML/level.5/tensile
python3 generate_tensile_mesh.py             # writes tensile.geom
../../../build/bin/tahoe -f tensile.xml      # ~30 s serial
python3 compare_to_hardening.py              # writes CSV + 2 PNGs
```

## Result — radial-return tracks the input spline within ~2 %

Probe IP at the centre of the bar (x = y = 0, z = 2.5 mm).  The bar is
in clean uniaxial tension there: `σ_xx`, `σ_yy ≈ 0` (≤ 10⁻⁵ MPa
throughout), so `σ_VM = σ_zz` and the current yield surface is
`σ_y(α) = σ_VM`.

| t | α | σ_zz [MPa] | σ_VM [MPa] | nearest knot (α, σ_y) |
|---:|----:|-----------:|-----------:|-----------------------|
| 0.04 | 0.003 | 288.6 | 288.7 | (0.005, 320)  — pre-knot |
| 0.20 | 0.018 | 410.9 | 411.3 | (0.020, 420) |
| 0.40 | 0.037 | 472.7 | 473.1 | (0.050, 500) — between (0.020, 420)·(0.050, 500) |
| 0.60 | 0.056 | 509.3 | 509.8 | (0.050, 500) |
| 0.80 | 0.074 | 528.7 | 529.3 | (0.080, 535) |
| 1.00 | 0.093 | 546.0 | 546.5 | between (0.080, 535)·(0.110, 560) — linear interp = 545.8 |

The simulation tracks the cubic spline within **~2 %** across the full
ε_p = 0 → 0.09 range — well within the radial-return tolerance.  The
slightly-low values mid-knot are the spline curvature (cubic dips below
linear interpolation between widely-spaced knots, by design).

![hardening](tensile_hardening.png)
![P-delta](tensile_pdelta.png)

## Why we don't pull to ε_p = 0.5 (the last spline knot)

The original plan was to pull to ε ≈ 0.5 to validate the full spline.
At about α = 0.10 the spline saturates: dσ_y/dε_p drops from ~800 MPa
to ~200 MPa over the 0.15–0.50 range.  With prescribed-displacement
control, that saturation makes the global Newton tangent become near-
singular in the axial direction (the change in axial force per unit
axial displacement collapses).  Line search reduces to step factor ≈ 0
and the iteration stalls.  Reaching the saturated portion of the curve
needs either arc-length control or a load-controlled formulation —
deferred as future work.  The 0 → 0.1 range tested here covers the
steepest part of the hardening response, which is where verification
matters most.

## See also
- `level.5/brinell/` — same material under contact + plasticity
- Issue #49 — tracking
- `tahoe/src/elements/continuum/solid/materials/plasticity_J2/J2_C0HardeningT.{h,cpp}` — the C¹-function-driven yield-stress interface used by `<Simo_J2>`
