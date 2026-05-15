# Level 3 — MPI distributed-memory parallel

Same kernels as levels 0–2, exercised under MPI domain decomposition with
4 ranks.  Requires `-DTAHOE_MPI=ON` and system OpenMPI (`libopenmpi-dev`);
do *not* use a conda MPI wrapper — see the top-level [`ReadMe`](../ReadMe).

Run:
```bash
cd benchmark_XML/level.3/parallel
/usr/bin/mpirun -np 4 ../../../build/bin/tahoe -f run.batch
# compare (serial) against per-rank reference output:
for x in *.xml; do ../../../build/bin/compare -f "$x"; done
```
The repo-root `run_benchmarks.sh level.3` automates both steps.

## Subdirectories

| File / dir | Theme |
|------------|-------|
| `parallel/elastostatic*.xml` | Static large-strain solid, three solver variants: profile (default), Newton + line search (`N-LS`), preconditioned CG (`PCG`), PCG with K-factor recompute (`PCG.K`). |
| `parallel/explicit.dyn{.restart}?.xml` | Explicit central-difference dynamics, with and without restart-file write/read. |
| `parallel/2D.particle/`, `3D.particle/` | Particle MD across ranks (cell-list communication, ghost atoms). |
| `parallel/2D.particle.spatial/`, `3D.particle.spatial/` | Particle MD with spatially-decomposed neighbour search. |
| `parallel/split_io/` | Multi-file output verification across ranks. |
| `parallel/geometry/` | Shared `.geom` files (decomposition output `*.nN.pK.geom*` is git-ignored). |

## Status (May 2026)

22 / 22 PASS, 0 FAIL.  Level.3 is currently clean.
