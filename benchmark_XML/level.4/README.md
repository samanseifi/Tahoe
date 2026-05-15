# Level 4 — Solver performance

Performance benchmarks for comparing the bundled sparse direct solvers.
Unlike levels 0–3 (correctness), level.4 problems are sized to stress the
linear solver with non-trivial DOF counts and multi-physics coupling.

## Subdirectories

| Directory | Theme |
|-----------|-------|
| `liquid_inclusion/` | Dielectric elastomer with a circular hole — 4704 Q4, 18 832 DOFs (displacement + electric potential).  The headline solver-comparison case, detailed below. |
| `single_element/` | Single-Hex8 sanity tests (monolithic vs staggered electro-mechanics; static and dynamic). |
| `dielectric_elastomer_damage/` | Three-way and two-way coupled staggered runs for dielectric elastomer + phase-field damage. |
| `surface_tension/` | 3D surface-tension Q1P0 patch test. |
| `geometry/` | Shared `.geom` meshes referenced by the cases above (partitioned `*.nN.pK.geom*` files are git-ignored). |

---

## Headline benchmark: dielectric elastomer with liquid inclusion (`liquid_inclusion/`)

A 2D coupled electromechanical problem using the `dielectric_elastomer_Q1P0Elastocapillary` element (mixed Q1P0 pressure + electric scalar potential field). The mesh is a square domain with a circular hole (`ah025.geom`).

| Property | Value |
|---|---|
| Mesh | `geometry/hole/ah025/ah025.geom` |
| Elements | 4704 Q4 quadrilateral elements |
| DOFs | **18,832** (displacement + electric potential) |
| Steps | 5 (benchmark subset of the 1470-step full run) |
| Newton iters / step | ~6 (convergence to relative tol 1e-9) |
| Physics | Large-strain electromechanics, incompressible pressure stabilisation |

Full problem: `voltage1_mumps.xml` (1470 steps). Benchmark variants use 5 steps and suppress output.

---

## Solver Benchmark Results

Measured on **WSL2 (Ubuntu 22.04), 12 logical cores** (Intel), March 2026.
Build: `cmake -B build -DTAHOE_MUMPS=ON -DTAHOE_MPI=ON -DTAHOE_SUPERLU=ON -DTAHOE_SPOOLES_MT=ON -DTAHOE_METIS=ON`.
Each entry is the **median of 3 trials** (wall-clock time, seconds).

### Serial / Single-Node Solvers

| Solver | XML element | Build flag | Median (s) | vs SPOOLES |
|--------|------------|------------|-----------|------------|
| **MUMPS** | `<MUMPS_matrix/>` | `TAHOE_MUMPS=ON` | **7.47** | **1.80× faster** |
| SPOOLES | `<SPOOLES_matrix/>` | default | 13.44 | baseline |
| SPOOLES-MT 4 threads | `<SPOOLES_MT_matrix num_threads="4"/>` | `TAHOE_SPOOLES_MT=ON` | 21.23 | 1.58× slower |
| SuperLU | `<SuperLU_matrix/>` | `TAHOE_SUPERLU=ON` | 26.08 | 1.94× slower |
| SPOOLES-MT 8 threads | `<SPOOLES_MT_matrix num_threads="8"/>` | `TAHOE_SPOOLES_MT=ON` | 28.89 | 2.15× slower |

> SPOOLES-MT and SuperLU are slower than serial SPOOLES on this problem due to WSL2 thread-scheduling overhead and the moderate problem size (18k DOFs). Multithreaded solvers scale better on bare-metal Linux with larger problems.

### MPI Distributed Solvers — Built-in Partitioner

Run with `/usr/bin/mpirun -np N ./build/bin/tahoe -f input.xml -decomp_method -0`.
Geometry pre-partitioned using the built-in recursive bisection partitioner.

| Solver | Ranks | XML element | Median (s) | vs serial MUMPS |
|--------|-------|------------|-----------|----------------|
| **MUMPS-MPI** | **2** | `<MUMPS_MPI_matrix/>` | **5.16** | **1.45× faster** |
| MUMPS-MPI | 4 | `<MUMPS_MPI_matrix/>` | 6.38 | 1.17× faster |
| MUMPS-MPI | 6 | `<MUMPS_MPI_matrix/>` | 6.78 | 1.10× faster |
| SPOOLES-MPI | 2 | `<SPOOLES_matrix/>` | 9.67 | 1.28× slower |
| SPOOLES-MPI | 4 | `<SPOOLES_matrix/>` | 9.28 | 1.24× slower |
| SPOOLES-MPI | 6 | `<SPOOLES_matrix/>` | 7.80 | 1.04× slower |

### MPI Distributed Solvers — METIS 5 Partitioner (`-DTAHOE_METIS=ON`)

Same runs with geometry re-partitioned by METIS 5 (`METIS_PartGraphKway`, edge-cut objective).
METIS is used automatically when built with `TAHOE_METIS=ON`; disable at runtime with `-no_metis`.

| Solver | Ranks | XML element | Median (s) | vs built-in | vs serial MUMPS |
|--------|-------|------------|-----------|------------|----------------|
| **MUMPS-MPI** | **2** | `<MUMPS_MPI_matrix/>` | **4.36** | **+16% faster** | **1.71× faster** |
| MUMPS-MPI | 4 | `<MUMPS_MPI_matrix/>` | 5.65 | +11% faster | 1.32× faster |
| MUMPS-MPI | 6 | `<MUMPS_MPI_matrix/>` | 6.19 | +9% faster | 1.21× faster |
| SPOOLES-MPI | 2 | `<SPOOLES_matrix/>` | 15.07 | — (high WSL2 variance) | 2.02× slower |
| SPOOLES-MPI | 4 | `<SPOOLES_matrix/>` | 12.47 | — (high WSL2 variance) | 1.67× slower |
| SPOOLES-MPI | 6 | `<SPOOLES_matrix/>` | 13.13 | — (high WSL2 variance) | 1.76× slower |

> METIS produces better-balanced partitions with fewer boundary DOFs, which directly reduces
> the inter-rank communication volume in MUMPS-MPI (distributed assembled solve). The improvement
> is most pronounced at np=2 (16%). SPOOLES-MPI timing shows high variance on WSL2 due to
> OS scheduler interference; results on bare-metal Linux are expected to be more consistent
> and also benefit from METIS partitioning.

---

## Overall Winner (with METIS)

```
MUMPS-MPI np=2 + METIS   4.4 s   ← fastest overall
MUMPS-MPI np=2 built-in  5.2 s
MUMPS serial             7.5 s   ← fastest single-process
SPOOLES-MPI np=6         7.8 s   (built-in; METIS variance too high on WSL2)
SPOOLES serial          13.4 s   ← default solver
```

For single-node production runs: use `<MUMPS_matrix/>` (serial MUMPS, no `mpirun` needed).
For MPI runs: use `<MUMPS_MPI_matrix/>` with 2 ranks and `-DTAHOE_METIS=ON`.

---

## Benchmark Input Files

| File | Solver | Notes |
|------|--------|-------|
| `voltage1_bench_spooles.xml` | SPOOLES (serial/MPI) | Use with `mpirun` for SPOOLES-MPI |
| `voltage1_bench_mumps.xml` | MUMPS serial | Single process |
| `voltage1_bench_mumps_mpi.xml` | MUMPS MPI | Use with `mpirun -decomp_method -0` |
| `voltage1_bench_superlu.xml` | SuperLU | Requires `TAHOE_SUPERLU=ON` |
| `voltage1_bench_spooles_mt4.xml` | SPOOLES-MT 4t | Requires `TAHOE_SPOOLES_MT=ON` |
| `voltage1_bench_spooles_mt8.xml` | SPOOLES-MT 8t | Requires `TAHOE_SPOOLES_MT=ON` |
| `voltage1_mumps.xml` | MUMPS serial | Full 1470-step production run |

## How to Run

```bash
TAHOE=/path/to/build/bin/tahoe
cd benchmark_XML/level.4/liquid_inclusion

# Serial solvers
$TAHOE -f voltage1_bench_mumps.xml
$TAHOE -f voltage1_bench_spooles.xml
$TAHOE -f voltage1_bench_superlu.xml

# MPI solvers — built-in partitioner (requires TAHOE_MPI=ON)
/usr/bin/mpirun -np 2 $TAHOE -f voltage1_bench_mumps_mpi.xml -decomp_method -0
/usr/bin/mpirun -np 4 $TAHOE -f voltage1_bench_spooles.xml   -decomp_method -0

# MPI solvers — METIS partitioner (requires TAHOE_MPI=ON + TAHOE_METIS=ON)
# Delete existing .n*.p*.geom files to force METIS re-decomposition
/usr/bin/mpirun -np 2 $TAHOE -f voltage1_bench_mumps_mpi.xml -decomp_method -0
/usr/bin/mpirun -np 4 $TAHOE -f voltage1_bench_mumps_mpi.xml -decomp_method -0
/usr/bin/mpirun -np 6 $TAHOE -f voltage1_bench_mumps_mpi.xml -decomp_method -0
```
