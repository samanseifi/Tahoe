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

## Status (May 2026)

47 / 47 PASS, 0 FAIL.  Level.2 is currently clean.
