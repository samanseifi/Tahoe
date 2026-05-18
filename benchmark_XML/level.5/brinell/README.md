# Brinell indentation — rigid ball on elastoplastic block (#47)

Quasi-static implicit indentation of a near-rigid hemispherical indenter
into a J2-plastic block, with light Coulomb friction.  Same quarter-
symmetry topology as the Hertz benchmark (`level.5/hertz/`), but the
base block is elastoplastic and the indenter is ~10× stiffer, so plastic
flow concentrates in the block.

The benchmark exercises three capabilities together for the first time:

| Capability | First landed in | Used here as |
|------------|-----------------|--------------|
| `cubic_spline` isotropic hardening for J2 | original Tahoe | takes 6 (ε_p, σ_y) points → C¹ spline through them |
| Implicit Coulomb friction in `PenaltyContact3DT` | PR #45 (#40) | μ = 0.1 sliding friction at the contact pair |
| Newton + line search (`<nonlinear_solver_LS>`) | original Tahoe | required for robust convergence under changing contact set + plasticity |

## Files

| File | Notes |
|------|-------|
| [`generate_brinell_mesh.py`](generate_brinell_mesh.py) | Writes both mesh variants (smoke, fine). |
| [`brinell_smoke.xml`](brinell_smoke.xml) | **Smoke variant** — 2 304 Hex8, 10 steps to δ = 0.15 mm, converges serial in ~55 s.  This is what CI / autonomous loops should run. |
| [`brinell.xml`](brinell.xml) | Fine variant — 8 800 Hex8, 20 steps to δ = 0.30 mm.  Slower (~1 h serial); deeper into the fully-plastic regime, suitable for Tabor-relation validation. |
| `brinell_smoke.geom`, `brinell.geom` | Generated meshes (text Tahoe geom format). |

## Setup

Both variants share the same geometry shape and BCs:

```
Indenter (block 1): hemisphere R = 5 mm, quarter-symmetry,
                    near-rigid elastic (E = 2 000 000 MPa, ν = 0.3)
Block    (block 2): R_box × R_box × H_base = 2.5 × 2.5 × 3 mm,
                    J2 plastic, E = 200 000 MPa, ν = 0.3,
                    σ_y0 = 250 MPa
                    cubic-spline hardening through:
                        (ε_p, σ_y) = (0.000, 250), (0.005, 320),
                                     (0.020, 420), (0.050, 500),
                                     (0.150, 580), (0.500, 650)  [MPa]
```

Boundary conditions:
- indenter top (NS1): prescribed `u_z = -δ_max · schedule(t)`
- indenter and block symmetry on x=0 and y=0 faces
- block bottom (NS4): fully clamped
- contact pair: `<contact_3D_penalty>` with `μ = 0.1`, regularised slip
  scale `friction_epsilon_velocity = 1e-4`, `penalty_stiffness = 1e7`

Solver: `<nonlinear_solver_LS>` + SPOOLES, `max_iterations=40`,
`search_iterations=5`.

## Running

```bash
cd benchmark_XML/level.5/brinell
python3 generate_brinell_mesh.py        # writes brinell_smoke.geom + brinell.geom
../../../build/bin/tahoe -f brinell_smoke.xml   # ~55 s on serial SPOOLES
# or for the fine variant:
../../../build/bin/tahoe -f brinell.xml         # ~1 h
```

## Convergence (brinell_smoke.xml)

All 10 steps converge.  Step 10 (deepest load, δ = 0.15 mm, fully plastic):

```
init: 0 LS: 3 ... | 0: Rel error = 1.41e-02
                   | 1: Rel error = 1.60e-02
                   | 2: Rel error = 9.07e-04
                   | 3: Rel error = 6.76e-04
                   | 4: Rel error = 3.57e-04
                   | 5: Rel error = 5.93e-05
                   | 6: Rel error = 2.24e-06
                   | 7: Rel error = 1.21e-08    ← converged
```

Contact patch grows monotonically as plastic deformation accumulates:
48 active strikers at step 1 → 80 at step 9 (plastic impression spreading).

## Expected physics (validated by `brinell.xml` on the fine mesh)

1. **Elastic phase** (δ ≲ 0.02 mm) — `P-δ` matches Hertz `P = (4/3) E* √R δ^{3/2}`.
2. **Yield onset** — maximum von Mises beneath the indenter reaches σ_y0 = 250 MPa
   at indentation depth `δ_y ≈ 0.012 mm` (Hertz analytical `p₀,y = 1.6 σ_y0`).
3. **Plastic phase** — contact patch grows faster than Hertz; mean pressure
   `p_m = P / πa²` flattens.
4. **Tabor relation** — at full load (δ = 0.3 mm in the fine variant),
   `HB ≈ p_m ≈ 2.8 σ_y` (Tabor 1951 for fully-plastic indentation).

## See also
- Issue #47 — Brinell benchmark tracking
- `level.5/hertz/` — same topology, elastic-only Hertz validation
- `level.5/implicit_friction/sliding_cubes.xml` — bare Coulomb friction test (#40)
- `level.5/tet_classic/tet4_hyperelastic_anp.xml` — implicit Newton + ANP-Tet4 (#29)

## Limitations

- The "rigid" indenter is actually a stiff elastic body (E = 2 × 10⁶ MPa).
  Its compliance contributes ~0.5 % of total displacement.  A true rigid-body
  element would be cleaner — tracked as #30.
- The fine variant (`brinell.xml`) takes ~1 h serial on profile_matrix
  with line-search Newton.  Each step is ~3 min, dominated by SPOOLES
  factorisation of the 29 022-DOF tangent.  For repeated runs, MUMPS or
  SuperLU would speed this up; not changed here so the benchmark is
  reproducible on a default build with no optional flags.
- Unloading / residual depth is not in this benchmark; would require an
  additional load schedule and a small code-side change to record
  per-step results without rewriting from scratch.
