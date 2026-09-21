# Level 0 — Core capabilities

Foundational regression tests covering every primary kernel: springs, solid
mechanics, diffusion, particles, meshfree, time integrators, I/O.  All code
changes should pass these before committing.

Run from the level root:
```bash
cd benchmark_XML/level.0
printf "run.batch\nquit\n" | ../../build/bin/tahoe    # run the suite
printf "run.batch\nquit\n" | ../../build/bin/compare  # compare to references
```

## Subdirectories

| Directory | Theme |
|-----------|-------|
| `2D.elastostatic/`, `3D.elastostatic/` | Static linear and large-strain solid mechanics; patch and CSE tests |
| `2D.elastodynamic/`, `3D.elastodynamic/` | Dynamic solid mechanics, central-difference and Newmark integrators |
| `2D.spring/`, `3D.spring/` | Truss and spring elements (smallest models in the suite) |
| `2D.particle/`, `3D.particle/` | Particle methods (Lennard-Jones, harmonic, EAM smoke tests) |
| `adhesion/` | Surface-to-surface adhesion (needs `ADHESION_ELEMENT` — enabled by default since #37) |
| `axisymmetric/` | Axisymmetric solid mechanics |
| `bridging/` | Atomistic-continuum coupling (needs `BRIDGING_ELEMENT` — currently off; runtime crash tracked in #37) |
| `diffusion/` | Heat / mass transport, linear and hyperbolic |
| `enhanced.strain/` | Enhanced-strain element formulations |
| `inputoutput/` | I/O round-trip tests |
| `integrator/` | Time-integrator verification (verlet, Newmark, HHT-α, …) |
| `matrix_check/` | Linear-solver sanity checks against analytical inverses |
| `meshfree/` | EFG / reproducing-kernel meshfree methods |
| `meshfree_kl_shell/` | Meshfree RKPM Kirchhoff–Love shell (#59): elastic obstacle course (Scordelis–Lo, hemisphere, pinched cylinder, #66) and essential-BC enforcement on the non-interpolatory field — collocation KBC (#70) and `penalty_displacement_meshfree` on the explicit and implicit paths (#73). |
| `phase_field/` | Phase-field fracture smoke tests |
| `geometry/` | Shared `.geom` mesh files referenced by tests in sibling directories |

## Status (May 2026, default build)

September 2026: `meshfree_kl_shell/` grew to 6 decks (#73). `scordelis_lo.xml` currently FAILS its
reference (rel. error 18 %): the reference was locked under #66 before the element's kernel default
moved to the paper's cubic B-spline with dilation 1.0 (#64), which gives −0.194 on this coarse
225-node cloud versus −0.27 with dilation 1.4 (see `docs/meshfree_kl_shell/history/`). Deliberately
left failing until the kernel/dilation choice for the coarse deck is settled; the other five pass.

172 / 186 PASS, 14 FAIL — all 14 remaining failures are real Tahoe code
regressions tracked under issue #37 (9 bridging runtime crashes, plus
`CSE.2.xml`, `adhesion.{2,3}.xml`, `inputoutput/square.xml`).  See
[`../ReadMe`](../ReadMe) for the full status table.
