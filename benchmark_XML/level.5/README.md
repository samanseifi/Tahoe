# Level 5 — Tahoe explicit-solver benchmarks

Benchmarks introduced alongside the modernised explicit-solver track
(`ExplicitElementT`, MVSIZ-batched).  Each subdirectory has its own
focus; tests reference these for performance baselines and physics
verification.

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

### Contact (issues #19, #26)

| File | Notes |
|------|-------|
| `vectorized_cubes_contact.xml` | Static compression, prescribed displacement (#19) |
| `vectorized_cubes_impact.xml` | Dynamic impact, initial velocity (#19) — bouncing |
| `vectorized_cubes_nofric.xml` | Sliding cube, frictionless (μ = 0) — reference |
| `vectorized_cubes_friction.xml` | Sliding cube, μ = 0.3 — Coulomb friction retards motion 10 % at t = 0.5 µs (#26) |

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

## tet_classic/

Classic-Tahoe (UpdatedLagrangianT) implicit Tet4 path — companion to
the explicit benchmarks above.

| File | Element | Material | Notes |
|------|---------|----------|-------|
| `tet4_hyperelastic.xml` | `<updated_lagrangian>` Tet4 | Neo-Hookean | Quasi-static stretch, classic Tet4 |
| `tet4_hyperelastic_anp.xml` | `<bonet_tet>` Tet4 | Neo-Hookean | Same with ANP F-bar (issue #28) |
| `compress_tet4.xml` | Plain Tet4 | Near-incompressible NH | Tet4 locks (κ = 1000, μ = 1) |
| `compress_bonet_tet.xml` | BonetTet | Same NH | ANP residual / inherited tangent inconsistency → issue #29 |

> **Known limitation (BonetTetT, classic implicit only)**: the ANP
> residual uses F̄ but the analytical tangent inherited from
> `UpdatedLagrangianT` is computed from F.  Newton stalls under
> near-incompressibility — fix is numerical tangent (forward
> differences of `FormKd`), tracked in issue #29.

## Coverage matrix

| Issue | Headline benchmark / test |
|-------|---------------------------|
| #16 J2 plasticity batch | `vectorized_hex_small_j2.xml`, `tests/materials/test_ExplJ2Plasticity.cpp` |
| #18 Taylor bar | `taylor_bar/taylor_ci.xml`, `tests/benchmarks/test_TaylorBar.cpp` |
| #19 Contact + explicit | `vectorized_cubes_contact.xml`, `vectorized_cubes_impact.xml` |
| #26 Coulomb friction | `vectorized_cubes_friction.xml`, `tests/benchmarks/test_FrictionContact.cpp` |
| #27 Tet4 kernel + Ji-fix | `vectorized_tet_*.xml`, `tests/kernels/test_Tet4Kernel.cpp` |
| #28 ANP-Tet4 (ELFORM=13) | `vectorized_tet_anp_small.xml`, `tests/kernels/test_ANPHelperT.cpp`, `tet_classic/*` |
