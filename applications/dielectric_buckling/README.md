# Dielectric-elastomer plate buckling (3D)

A slender 3D dielectric-elastomer plate with clamped ends and side rollers,
loaded by a transverse electric field. The Maxwell stress squeezes the plate
through its thickness; Poisson coupling pushes it laterally; clamped ends
convert that into compressive in-plane stress, which ultimately buckles the
plate out of plane. A small mode-1 sinusoidal imperfection is baked into the
mesh to seed the instability.

We compare two time-integration strategies on the same staggered Maxwell
coupling: **explicit central-difference** vs **implicit nonlinear-HHT**. The
quasi-static buckled shape should be the same; the explicit run takes 5–10×
more steps but each step is cheap.

## Geometry

| | length (X) | width (Y) | thickness (Z) |
|---|---|---|---|
| size | 20 | 2 | 1 |
| elements | 80 | 8 | 4 |

- Slenderness L/t = 20 (clamped-clamped Euler regime).
- Imperfection: `z += 0.01 * sin(pi * x / Lx)` on every node — vanishes at the
  clamps, peaks at midspan.

Boundary conditions:

| Node set | Location | Mechanics | Electrical |
|---|---|---|---|
| 1 | x = 0 | clamped (fix all 3 DOFs) | — |
| 2 | x = Lx | clamped (fix all 3 DOFs) | — |
| 3 | y = 0 | roller (fix D_Y) | — |
| 4 | y = Ly | roller (fix D_Y) | — |
| 5 | z = 0 | — | Psi = 0 |
| 6 | z = Lz | — | Psi = V(t) |

Voltage ramps linearly from 0 to 1.5 over T = 2.0, then holds.

## Material

Neo-Hookean dielectric elastomer (dimensionless), mildly compressible to
avoid volumetric locking with the standard hex8:

| | | |
|---|---|---|
| mu | 1.0 | shear modulus |
| lambda | 100 | first Lamé parameter (bulk modulus kappa = 100.67 ≈ 100·mu) |
| epsilon | 1.0 | dielectric permittivity |
| density | 1.0 | |

## Two simulations

Both XMLs share the same mesh, BCs, voltage schedule, and Maxwell-stress
coupling — only the mechanical time integrator differs.

### `staggered_implicit.xml`
- **Pass 1 — implicit static electrostatics**: `<diffusion>` element with
  `linear_dielectric_material` solves `div(eps grad Psi) = 0` via SPOOLES.
- **Pass 2 — implicit dynamic mechanics**: `<updated_lagrangian_Q1P0>` with
  `nonlinear_HHT` integrator and Newton-Raphson via SPOOLES. The element
  auto-detects the `electric_scalar_potential` field and adds the Maxwell
  stress `sigma_E = eps (E ⊗ E - 1/2 |E|^2 I)` to its residual.

`dt = 0.01`, 200 steps. Unconditionally stable, ~50 s wall time.

### `staggered_explicit.xml`
- **Pass 1**: same as above.
- **Pass 2 — explicit dynamic mechanics**: same element with the
  `central_difference` integrator. Lumped mass, no matrix solve — each step
  is a diagonal divide of the residual.

`dt = 0.005` (CFL on `c_p = sqrt((lambda+2*mu)/rho) ≈ 10.1`, `h_min = 0.25`),
400 steps. ~11 s wall time.

## Why not a monolithic 4-DOF/node reference?

The plan was to use one of the bundled fully-coupled 3D dielectric elements
as the "ground truth", but both 3D variants are unusable in this build:

- **`dielectric_elastomer_Q1P0` (3D Q1P0)** — tangent goes NaN on the very
  first iteration even on a single hex. The material registers
  `Young_Modulus` and `Poisson` parameters that are never used in the
  stress/tangent code, and there is **no benchmark XML anywhere in the tree
  that exercises it**. The 3D Q1P0 path appears unfinished.
- **`dielectric_elastomer` (3D, non-Q1P0)** — works on a single hex; on a
  multi-element mesh both `SPOOLES_matrix` (`double free or corruption`) and
  `MUMPS_matrix` (`malloc(): invalid size`) abort during assembly of the
  asymmetric coupled 4-DOF/node tangent. `profile_matrix` works but is
  unusably slow at 14k DOFs.

So the only validated 3D dielectric path is the staggered one
(`<diffusion>` + `<updated_lagrangian_Q1P0>` with electric-field coupling
through `SimoQ1P0`'s Maxwell-stress hook), which is what both XMLs use.

## Running

```bash
# 1. Generate the mesh
python3 generate_plate_mesh.py

# 2. Run either case
../../build/bin/tahoe -f staggered_implicit.xml
../../build/bin/tahoe -f staggered_explicit.xml
```

Outputs are ExodusII (`*.io1.exo`); open in ParaView and look at `D_Z` along
the midline (`y = Ly/2, z = Lz/2`) to compare the buckled shapes.

## Validation status (V=1.5, T=2.0 baseline)

| field | implicit range | explicit range | max abs diff | relative |
|---|---|---|---|---|
| D_X | ±2.41e−3 | ±2.41e−3 | 5.2e−6 | 0.2% |
| D_Z | [−1.25e−2, +1.07e−2] | [−1.24e−2, +1.06e−2] | 1.9e−5 | 0.2% |
| s11 | [−2.34, −2.26] | [−2.34, −2.26] | 6.3e−4 | <0.1% |

Maxwell thickness squeeze: -8.5620e−3 vs -8.5619e−3 (5-digit agreement).
Wall time: implicit 141 s, explicit 11 s.

**The plate is *not* fully buckled at V=1.5/T=2.0.** Compressive stress σ11 ≈ -2.3
is far past the linearized Euler critical (~0.025), but without damping the
dynamic explicit/implicit-HHT trajectories don't relax to the post-buckled
quasi-static state in the available time. The imperfection does grow
exponentially (`5.8e−5 → 1.6e−3 → 2.5e−3` over t = 1.5→2.8 in a V=2.5/T=4
trial), so longer simulation + damping or an arc-length static solver is
needed to reach the visible post-buckled state. The baseline runs are good
**integrator validation**, not a finished buckling study.

### What was tried, and why a clean buckle is still open

- **Higher V (2.5–3.0) + bigger imperfection (5%) explicit dynamic**: ran
  for ~480 of 1600 steps, then hit `zero or negative jacobian` — element
  inversion at the clamped corners (e.g. node at `(0.75, *, 1.006)` and
  `(19.25, *, 1.006)`) where the thickness squeeze concentrates against
  the rigid clamp.
- **Quasi-static `integrator="static"` Newton-Raphson, 50 increments**:
  same corner inversion at step 4 (V≈0.24).
- **Quasi-static, 200 increments, 1% imperfection**: corner inversion
  goes away, but now the *staggered* (Gauss-Seidel) coupling between the
  electrostatic and mechanical solves stagnates at step 34 (V≈0.51) —
  each individual phase converges tightly, but the outer staggered loop
  oscillates between them because the mech↔elec feedback is too strong
  at this load. Classic limitation of explicit-implicit / fixed-point
  staggering; would need under-relaxation, Aitken acceleration, or a
  monolithic solve to break the deadlock.
- **Implication**: getting a clean visible buckle out of this geometry
  needs one of:
   1. **Less aggressive clamping** at the end faces (e.g. only fix `D_X`
      on the end face, leaving `D_Y, D_Z` free) — sidesteps the corner
      squeeze.
   2. **A working monolithic 3D coupled element** — both
      `dielectric_elastomer_Q1P0` (NaN) and `dielectric_elastomer`
      (heap corruption with SPOOLES/MUMPS) are broken in this build.
   3. **An arc-length / branch-following continuation solver** — Tahoe
      doesn't appear to ship one for finite-strain coupled problems.

## Tweaking

- **Higher voltage / stronger buckle**: bump the `value="1.5"` on the top
  electrode `kinematic_BC` in both XMLs.
- **Slenderness**: rerun the mesh script with different `--Lx`. Critical
  voltage scales roughly as `1/L`.
- **Mesh refinement**: `--nx 160 --ny 16 --nz 8` doubles linear resolution
  (8× elements). Remember to halve `time_step` in the explicit case to keep
  the CFL condition.
- **Imperfection size**: `--imperfection 0.05` (5% of thickness) gives a
  stronger seed and reaches an obvious post-buckled state faster — useful
  for visual sanity checks before running with a tiny imperfection.
