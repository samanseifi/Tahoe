# Level 2 — Extended capabilities

Larger / longer-running problems than levels 0–1, exercising whole feature
combinations rather than single kernels.  No optional modules required.

Run as the other levels:
```bash
cd benchmark_XML/level.2
printf "run.batch\nquit\n" | ../../build/bin/tahoe
printf "run.batch\nquit\n" | ../../build/bin/compare
```

## Subdirectories

| Directory | Theme |
|-----------|-------|
| `K.field/` | Crack-tip K-field driven elastic / plastic patch (linear, large-strain, contour-integral variants). |
| `angled_bc/` | Skew kinematic boundary conditions (BCs not aligned with global axes). |
| `contact_simple/` | Frictionless penalty contact between simple Q4 / Hex8 bodies — the canonical contact regression. |
| `conveyor/` | Periodic-domain conveyor / moving-mesh tests; exercises restart files and the `*.tracking` output. |
| `force.controller/` | Prescribed-force controllers (PD, integrator-based) driving displacement BCs to a target traction. |
| `thermostats/` | Particle-thermostat schemes (Nose-Hoover, Langevin) on Lennard-Jones systems. |
| `tied/` | Tied-node constraints between dissimilar meshes. |
| `torsion/` | Torsional loading on solid sections, large rotation. |
| `geometry/` | Shared `.geom` files used by tests in sibling directories. |
| `meshfree_kl_shell/` | Elasto-plastic paper cases of the meshfree Kirchhoff–Love shell (#59, #73), shortened and mass-scaled: necking of a cylindrical shell (Wang & Bazilevs §4.3, J2 with saturation hardening, thickness update) and the elasto-plastic pinched cylinder (§4.4, `penalty_displacement_meshfree` loading, plastic fold). Each runs in ~45 s single-threaded; the full-resolution validation lives in `docs/meshfree_kl_shell/`. |

## Status (May 2026)

47 / 47 PASS, 0 FAIL.  Level.2 is currently clean.

September 2026: `meshfree_kl_shell/` added (2 decks, references generated with the merged develop of 2026-09-21).
