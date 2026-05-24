# Level 5 — Modernised-solver benchmarks

Benchmarks introduced alongside the modernised solver tracks — the explicit
`ExplicitElementT` (MVSIZ-batched, 2026 modernisation) and the matching
implicit additions on the classic `UpdatedLagrangianT` path (BonetTet ANP,
implicit Coulomb friction).  Each subdirectory has its own focus; tests
reference these for performance baselines and physics verification.

## Subdirectories at a glance

| Directory | Theme |
|-----------|-------|
| `explicit_benchmark/` | The original speed-up suite (legacy ↔ vectorized), plus growing feature-specific smoke tests (J2 plasticity, ANP-Tet4, contact, friction, hourglass control). |
| `taylor_bar/` | Classic Taylor-bar impact, soft-metal demo (issue #18). |
| `hertz/` | Hertz contact validation — Hex8-only and mixed Tet/Hex variants; explicit vs implicit cross-validated to within 0.5 %. |
| `tet_classic/` | Implicit Tet4 path: plain `<updated_lagrangian>` and ANP `<bonet_tet>` variants (#27, #28, #29). |
| `implicit_friction/` | Implicit Coulomb friction smoke test (#40) — residual (PR #41) + LHS tangent (PR #45) both landed; converges quadratically. |
| `brinell/` | Brinell-style indentation — rigid ball on J2-plastic block, light Coulomb friction (#47).  First benchmark combining cubic-spline hardening + implicit friction + Newton-LS. |
| `tensile/` | Uniaxial tension on a J2-plastic bar (#49) — verifies that `<Simo_J2>` with `cubic_spline` hardening reproduces the input σ_y(ε_p) data to ~2 % under monotonic loading. |


## explicit_benchmark/

The original speed-up suite (legacy vs MVSIZ vectorized) plus growing
families of feature-specific smoke tests.

### Performance comparison (Q4 / Hex8, the original measurements)

| Pair | Mesh | Purpose |
|------|------|---------|
| `legacy_small.xml` ↔ `vectorized_small.xml` | 200 Q4 | 2D speed baseline |
| `legacy_medium.xml` ↔ `vectorized_medium.xml` | 5 000 Q4 | 2D mid-size |
| `legacy_large.xml` ↔ `vectorized_large.xml` | 20 000 Q4 | 2D large |
| `legacy_hex_small.xml` ↔ `vectorized_hex_small.xml` | 250 Hex8 | 3D speed baseline |
| `legacy_hex_medium.xml` ↔ `vectorized_hex_medium.xml` | 2 000 Hex8 | 3D mid-size |
| `legacy_hex_large.xml` ↔ `vectorized_hex_large.xml` | 8 000 Hex8 | 3D large |

Run with `bash run_benchmark.sh` (2D) or `bash run_3d_benchmark.sh`
(3D); reference timings live in [PERFORMANCE_ROADMAP.md](explicit_benchmark/PERFORMANCE_ROADMAP.md).

### Reduced-integration / hourglass / mass-scaling

| File | Element | Notes |
|------|---------|-------|
| `legacy_large_1ip.xml` | Q4 1-IP | Legacy reduced integration |
| `vectorized_large_1ip.xml` | Q4 1-IP | Vectorized reduced integration |
| `vectorized_large_1ip_hg.xml` | Q4 1-IP + viscous hourglass | Hourglass control demo |
| `vectorized_large_ms.xml` | Q4 + mass scaling | Adaptive mass-scaling demo |
| `vectorized_hex_*_1ip.xml` | Hex8 1-IP variants |  |

### J2 plasticity smoke (issue #16)

| File | Notes |
|------|-------|
| `vectorized_hex_small_j2.xml` | 250 Hex8 + J2 plasticity, soft-metal yield, end-to-end smoke |

### Tet4 element family (issue #27)

Increasing-complexity sequence used while bringing the Tet4 kernel up.

| File | Mesh | Purpose |
|------|------|---------|
| `vectorized_tet_one.xml` | 1 tet | Single-element rigid-body translation check |
| `vectorized_tet_two.xml` | 2 tets sharing a face | Shared-node assembly check |
| `vectorized_tet_hex.xml` | 5 tets (one hex split) | Element-group integrity |
| `vectorized_tet_minimal.xml` | 270 tets, free body | Mass / IC verification |
| `vectorized_tet_small.xml` | 270 tets, plain Tet4 + J2 + impact | Locked-volumetric reference |

### ANP-Tet4 (issue #28, LS-DYNA ELFORM=13 equivalent)

| File | Notes |
|------|-------|
| `vectorized_tet_anp_small.xml` | Same mesh as `vectorized_tet_small.xml` with `<anp_tet4 enabled="true"/>` — F-bar averaging produces ~35 % more compliance, exactly the locking relief expected from ELFORM=13 |

### Contact (issues #19, #26, #31, #32, #33)

| File | Notes |
|------|-------|
| `vectorized_cubes_contact.xml` | Static compression, prescribed displacement (#19) |
| `vectorized_cubes_impact.xml` | Dynamic impact, initial velocity (#19) — bouncing |
| `vectorized_cubes_nofric.xml` | Sliding cube, frictionless (μ = 0) — reference |
| `vectorized_cubes_friction.xml` | Sliding cube, μ = 0.3 — Coulomb friction retards motion 10 % at t = 0.5 µs (#26) |

Performance + stability fixes shipped alongside these benchmarks:
- `viscous_damping` parameter on `<contact_3D_penalty>` for normal-direction
  damping in explicit runs — quenches contact-mode ringing (#31)
- OpenMP auto-tune threshold (`num_batches >= 4 * num_threads`) on
  `ExplicitElementT` — single-thread performance preserved on small
  meshes that otherwise pay full OMP overhead (#32)
- OpenMP-parallel `RHSDriver` for `PenaltyContact3DT` with thread-local
  workspaces and `omp critical` only around the global-RHS scatter (#33).
  See `tests/benchmarks/test_ContactPerf.cpp` for the integration tests
  that exercise all three fixes together.

### Mesh generators

| Script | Output |
|--------|--------|
| `generate_mesh.py` | 2D Q4 (`block_small/medium/large`) |
| `generate_3d_mesh.py` | 3D Hex8 (`hex_small/medium/large`) |
| `generate_tet_mesh.py` | 3D Tet4 via 5-tet hex split (`tet_small/medium`) |

## taylor_bar/

Classic Taylor-bar impact, soft-metal demo parameters (issue #18).

| File | Runtime | Purpose |
|------|---------|---------|
| `taylor_ci.xml` | ~10 s on CI | Smoke test, 1 µs simulation |
| `taylor_small.xml` | ~10 min | Full 60 µs validation run |

Reference numbers in the commit message of `#18`; full Wilkins-Guinan
match would need a cylindrical geometry (current mesh is prismatic
quarter-symmetry, captures KE → plastic-work conversion correctly but
not the Hugoniot pressure concentration).

## hertz/

Quarter-symmetry curved-bottom indenter on an elastic base, sized so
both bodies sit cleanly in the half-space limit Hertz assumes.
Exercises both Tahoe contact pathways (explicit + implicit) on the
same mesh.  Match against analytical Hertz — `a`, `F`, `p₀`, `p_mid` —
within 1-2 % on every fitted metric.  See [`hertz/README.md`](hertz/README.md).

| File | Purpose |
|------|---------|
| `hertz_explicit.xml` | `<explicit_solid>` + central-difference + viscous damping |
| `hertz_implicit.xml` | `<updated_lagrangian>` + Newton-Raphson |
| `generate_hertz_mesh.py` | tanh-graded Hex8 mesh, 16 128 elements |
| `compare_to_analytical.py` | post-processing — extracts `a`, `F`, `p₀`, fits Hertz, plots |
| `hertz_pressure.png` | side-by-side dimensional + normalised p(r) overlay |

Cross-validation: explicit and implicit agree with each other to
better than 0.5 % on every fitted metric — strong evidence that
Tahoe's two contact pathways compute the same physics.

A **Tet4 indenter on Hex8 base** variant lives alongside it:
`hertz_tet_hex_*.xml` + `generate_hertz_tet_hex_mesh.py` +
`compare_tet_hex_to_analytical.py`.  Same geometry, same material, but
the indenter is meshed via the standard 5-tet split (46 080 Tet4)
driven by `<bonet_tet>` (implicit) or `<explicit_solid>` + `<anp_tet4>`
(explicit).  Demonstrates mixed-element-type contact (Tet-tri ↔
Hex-quad, with the quads auto-split to triangles by
`Contact3DT::ConvertQuadToTri`) and eliminates the +15 % apex-pressure
spike of the all-hex case.  `merge_io.py` combines the per-element-
group io files into one Exodus for ParaView.

## tet_classic/

Classic-Tahoe (UpdatedLagrangianT) implicit Tet4 path — companion to
the explicit benchmarks above.

| File | Element | Material | Notes |
|------|---------|----------|-------|
| `tet4_hyperelastic.xml` | `<updated_lagrangian>` Tet4 | Neo-Hookean | Quasi-static stretch, classic Tet4 |
| `tet4_hyperelastic_anp.xml` | `<bonet_tet>` Tet4 | Neo-Hookean (κ=100, μ=40) | ANP F-bar with lagged-J̄ + FD tangent (#29) — quadratic Newton, ~1 iter/step |
| `compress_tet4.xml` | Plain Tet4 | Near-incompressible NH (κ=1000, μ=1) | Tet4 locks |
| `compress_bonet_tet.xml` | BonetTet | Same NH | Diverges (severe near-incompressibility, see limitation below) |

> **#29 status — partial improvement landed.** The lagged-J̄ scheme
> caches J̄ at step start and a per-element forward-FD `FormStiffness`
> builds the tangent from the F̄ residual at frozen J̄.  On moderately
> near-incompressible cases (κ/μ ~ 2.5, `tet4_hyperelastic_anp.xml`)
> Newton now converges quadratically in ~1 iter/step — a 15× iteration
> reduction vs the linear-rate baseline that needed ~16 iter/step.
>
> **Remaining limitation**: the per-element 12×12 FD tangent cannot
> express the cross-element ∂J̄/∂u coupling (J̄ is nodal-averaged, so a
> node perturbation changes J̄ in every neighbouring element).  At
> severe near-incompressibility (κ ≫ μ, e.g. `compress_bonet_tet.xml`
> with κ=1000, μ=1) the missing block kills the local volumetric
> stiffness and Newton diverges.  Bonet & Burton 1998 itself only
> formulates the *residual* — for explicit dynamics — and never derives
> an implicit tangent; a fully consistent tangent would follow
> Bonet/Marriott/Hassan 2001 or Caylak/Mahnken 2012, kept open as a
> future analytical-tangent extension of #29.

## implicit_friction/

Quasi-static counterpart to the explicit Coulomb friction benchmark
(`explicit_benchmark/vectorized_cubes_friction.xml`, #26).  Same two-cube
geometry driven by `<updated_lagrangian>` + Newton-Raphson with the
slip-history Coulomb branch in `PenaltyContact3DT` (#40).

| File | Notes |
|------|-------|
| `sliding_cubes.xml` | Top face compressed (u_z = −0.005) and dragged tangentially (u_x = 0 → 0.01 over 100 steps); μ = 0.02. |
| `two_cubes.geom` | Copy of `../explicit_benchmark/geometry/two_cubes_refined.geom`. |

**Status: converges quadratically.**  The residual + slip-history landed
in PR #41; the LHS friction tangent (per-pair 12×12 forward-FD on
`ImplicitFrictionRHS`) landed in PR #45.  Newton converges in 1–2 iters
per step; whole 100-step run completes in ~4 s serial.  See
[`implicit_friction/README.md`](implicit_friction/README.md).

## brinell/

Quasi-static implicit indentation of a near-rigid hemispherical indenter
into a J2-plastic block, with light Coulomb friction (#47).  Same quarter-
symmetry topology as `hertz/` (curved-bottom block on flat half-space),
but the base is elastoplastic and the indenter is ~10× stiffer so plastic
flow concentrates in the block.

| File | Notes |
|------|-------|
| [`brinell_smoke.xml`](brinell/brinell_smoke.xml) | Smoke variant: ~2.3k Hex8, 10 steps, δ = 0 → 0.15 mm, ~55 s serial.  This is the CI-friendly path. |
| [`brinell.xml`](brinell/brinell.xml) | Fine variant: 8 800 Hex8, 20 steps, δ = 0 → 0.30 mm, ~1 h serial.  Goes deep enough into the fully-plastic regime to validate the Tabor relation `HB ≈ 2.8 σ_y`. |
| [`generate_brinell_mesh.py`](brinell/generate_brinell_mesh.py) | Mesh generator (writes both variants). |

Material: steel-like J2 with `cubic_spline` hardening through six
`(ε_p, σ_y)` points (250 → 650 MPa).  First benchmark combining
spline-hardening plasticity with implicit Coulomb friction (#40) and
Newton-LS (`<nonlinear_solver_LS>`).

## Coverage matrix

| Issue | Headline benchmark / test |
|-------|---------------------------|
| #16 J2 plasticity batch | `vectorized_hex_small_j2.xml`, `tests/materials/test_ExplJ2Plasticity.cpp` |
| #18 Taylor bar | `taylor_bar/taylor_ci.xml`, `tests/benchmarks/test_TaylorBar.cpp` |
| #19 Contact + explicit | `vectorized_cubes_contact.xml`, `vectorized_cubes_impact.xml` |
| #26 Coulomb friction | `vectorized_cubes_friction.xml`, `tests/benchmarks/test_FrictionContact.cpp` |
| #27 Tet4 kernel + Ji-fix | `vectorized_tet_*.xml`, `tests/kernels/test_Tet4Kernel.cpp` |
| #28 ANP-Tet4 (ELFORM=13) | `vectorized_tet_anp_small.xml`, `tests/kernels/test_ANPHelperT.cpp`, `tet_classic/*` |
| #29 BonetTet lagged-J̄ + FD tangent | `tet_classic/tet4_hyperelastic_anp.xml` (quadratic on moderate κ) |
| #31/#32/#33 contact-stack perf | `tests/benchmarks/test_ContactPerf.cpp`, `vectorized_cubes_*.xml` (viscous damping, OMP threshold, parallel contact) |
| #40 implicit Coulomb friction | `implicit_friction/sliding_cubes.xml` (residual + LHS tangent — converges, quadratic Newton) |
| #47 Brinell indentation       | `brinell/brinell_smoke.xml` (~55 s) + `brinell/brinell.xml` (Tabor validation) |
| #49 J2 spline-hardening verification | `tensile/tensile.xml` — σ_VM(α) overlays input spline within ~2 % |
| Hertz contact validation | `hertz/hertz_explicit.xml`, `hertz/hertz_implicit.xml`, `hertz/compare_to_analytical.py` |
| Mixed Tet/Hex contact | `hertz/hertz_tet_hex_*.xml`, `hertz/compare_tet_hex_to_analytical.py` |
