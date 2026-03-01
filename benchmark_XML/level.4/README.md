# Level 4 Benchmarks — Solver Performance

Level 4 contains performance benchmarks for comparing sparse direct solvers. Unlike levels 0–3 (which test correctness), level 4 problems are designed to stress the linear solver with non-trivial DOF counts and multi-physics coupling.

---

## Problem: Dielectric Elastomer with Liquid Inclusion (`liquid_inclusion/`)

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
Build: `cmake -B build -DTAHOE_MUMPS=ON -DTAHOE_MPI=ON -DTAHOE_SUPERLU=ON -DTAHOE_SPOOLES_MT=ON`.
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

### MPI Distributed Solvers

Run with `/usr/bin/mpirun -np N ./build/bin/tahoe -f input.xml -decomp_method -0`.

| Solver | Ranks | XML element | Median (s) | vs serial MUMPS |
|--------|-------|------------|-----------|----------------|
| **MUMPS-MPI** | **2** | `<MUMPS_MPI_matrix/>` | **5.16** | **1.45× faster** |
| MUMPS-MPI | 4 | `<MUMPS_MPI_matrix/>` | 6.38 | 1.17× faster |
| MUMPS-MPI | 6 | `<MUMPS_MPI_matrix/>` | 6.78 | 1.10× faster |
| SPOOLES-MPI | 2 | `<SPOOLES_matrix/>` | 9.67 | 1.28× slower |
| SPOOLES-MPI | 4 | `<SPOOLES_matrix/>` | 9.28 | 1.24× slower |
| SPOOLES-MPI | 6 | `<SPOOLES_matrix/>` | 7.80 | 1.04× slower |

> MUMPS-MPI uses distributed assembled input (`icntl[17]=3`): each rank supplies its local COO triplets with global indices; MUMPS handles the distributed factorization. The RHS is gathered on rank 0 before the solve and broadcast back to all ranks. **2 ranks is the sweet spot** for this problem size — beyond that, MPI communication overhead limits further gains on WSL2.

---

## Overall Winner

```
MUMPS-MPI np=2   5.2 s   ← fastest
MUMPS serial     7.5 s   ← fastest single-process
SPOOLES-MPI np=6 7.8 s
SPOOLES serial  13.4 s   ← default solver
```

For single-node production runs: use `<MUMPS_matrix/>` (serial MUMPS, no `mpirun` needed).
For MPI runs (already using domain decomposition): use `<MUMPS_MPI_matrix/>` with 2–4 ranks.

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

# MPI solvers (requires TAHOE_MPI=ON build)
/usr/bin/mpirun -np 2 $TAHOE -f voltage1_bench_mumps_mpi.xml -decomp_method -0
/usr/bin/mpirun -np 4 $TAHOE -f voltage1_bench_spooles.xml   -decomp_method -0
```
