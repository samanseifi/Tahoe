# SPOOLES-MPI — Distributed-Memory Sparse Direct Solver

MPI parallel extension to SPOOLES 2.2. Distributes the LU factorisation
and solve across multiple MPI ranks, enabling large-scale parallel finite
element analysis on clusters. Sourced from the NguyenLabJHU/Tahoe-FEM fork,
which bundles the original `spoolesMPI` directory from the SPOOLES 2.2
distribution (Rice University).

## Requirements

| Dependency | Notes |
|------------|-------|
| `TAHOE_SPOOLES=ON` | Serial SPOOLES base library must be built |
| System OpenMPI | `sudo apt-get install libopenmpi-dev` on Debian/Ubuntu |
| `TAHOE_MPI=ON` | Enables both spoolesMPI and `__TAHOE_MPI__` in toolbox |

> **Important — system MPI vs conda MPI**: Tahoe's tahoe_spooles_mpi links
> against system OpenMPI (`libmpi.so.40` from `/lib/x86_64-linux-gnu/`).
> If conda is active, its `mpirun` uses MPICH, which is ABI-incompatible.
> Always use `/usr/bin/mpirun` (system OpenMPI) to launch the binary.
> See [Troubleshooting](#troubleshooting) below.

## Enabling

```bash
# Prerequisite (Debian/Ubuntu)
sudo apt-get install libopenmpi-dev

# Configure — uses --fresh to avoid stale conda MPI paths in cache
cmake -B build --fresh -DTAHOE_MPI=ON
cmake --build build -j$(nproc)
```

Verify the correct MPI was found:
```
-- Tahoe: MPI enabled  (/usr/bin/mpicxx)
-- Found MPI_C: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so (found version "4.1")
```

## Running

### Single input file
```bash
mpirun -np 4 ./build/bin/tahoe -f input.xml
```
Tahoe automatically decomposes the mesh across ranks, runs the analysis, and
writes per-rank output files (`input.io0.run`, `input.io1.run`, …) or a joined
output depending on I/O mode.

### Batch mode (for level.3 regression tests)
```bash
cd benchmark_XML/level.3/parallel
/usr/bin/mpirun -np 4 /path/to/build/bin/tahoe -f run.batch
```

Batch files use `@` as the first character (batch-mode marker). Option lines
(`-run`, `-join_io`, `-decomp_method -0`) set global state before each XML
file. Post-run comparison:
```bash
/path/to/build/bin/compare -f elastostatic.xml   # serial compare tool
```

### Decomposition options (set in batch or XML)

| Flag | Meaning |
|------|---------|
| `-decomp_method -0` | Sequential index decomposition (no METIS required) |
| `-decomp_method -1` | Spatial decomposition |
| `-no_metis` | Suppress graph-based decomposition (use with `-0` or `-1`) |
| `-join_io` | All ranks write joined output (default for top-level batch) |
| `-split_io` | Each rank writes `*.pN.ioM.run` files independently |
| `-join -N xml` | Merge N split outputs back into a single joined result |

## Source layout

```
spoolesMPI/src/
    MPI/
        MPI.factorMPI.c          distributed LU factorisation
        MPI.solveMPI.c           distributed triangular solve
        MPI.splitFrontMtx.c      frontal matrix distribution
        MPI.DenseMtx_gather.c    gather dense matrix data
        MPI.DenseMtx_scatterAdd.c scatter/add
        MPI.ETree_Bcast.c        broadcast elimination tree
        MPI.Graph_Bcast.c        broadcast graph
        MPI.IV_Bcast.c           broadcast integer vector
        MPI.IVL_Bcast.c          broadcast integer vector list
        MPI.IVL_alltoall.c       all-to-all for IVL
        MPI.IVLallgather.c       all-gather for IVL
        MPI.IVallgather.c        all-gather for IV
        MPI.MMM.c                distributed matrix-matrix multiply
        MPI.InpMtx_*.c           distributed input matrix operations
        MPI.FrontMtx_*.c         distributed frontal matrix operations
        … (25 files total)
        SPOOLESMPI.h             top-level include
    drivers/
        LU_MPI_driver.c          high-level entry point
        LU_MPI_driver.h          public API
        LU_MPI_driver_int.h      internal header
        LU_MPI_driver_init.c     initialise driver (read matrix, order, factor)
        LU_MPI_driver_factorize.c  numerical factorisation
        LU_MPI_driver_solve.c    forward/back solve
        LU_MPI_driver_free.c     free all structures
```

## Tahoe wrapper

```
tahoe/src/primitives/globalmatrix/SPOOLES/
    SPOOLESMatrixT_mpi.h / .cpp    distributed LU wrapper
```

`SolverT.cpp` selects `SPOOLESMatrixT_mpi` when both `__SPOOLES_MPI__` and
`__TAHOE_MPI__` are defined and the communicator size is greater than one.
For a single-rank MPI run it falls back to the serial `SPOOLESMatrixT`.

## Compile defines

| Define | Set by | Effect |
|--------|--------|--------|
| `__SPOOLES_MPI__` | `tahoe_spooles_mpi` (PUBLIC) | Activates MPI solver path in `SolverT.cpp` |
| `__TAHOE_MPI__` | `toolbox/CMakeLists.txt` (PUBLIC) | Activates `CommunicatorT`, `CommManagerT`, domain decomposition |

## MPI communication layer (`toolbox/src/parallel/`)

`CommunicatorT` is Tahoe's portable MPI abstraction. When `__TAHOE_MPI__` is
not defined it degrades to serial no-ops, so the same binary can be built
with or without MPI. Key classes:

| Class | File | Purpose |
|-------|------|---------|
| `CommunicatorT` | `CommunicatorT.h/.cpp` | Wraps `MPI_Comm`; broadcast, reduce, gather, scatter |
| `CommManagerT` | `CommManagerT.h/.cpp` | Manages inter-partition message schedules |
| `PartitionT` | `PartitionT.h/.cpp` | Mesh partitioning strategies |

## Benchmark results (level.3)

22/22 level.3 MPI regression tests pass with 4 MPI ranks using system
OpenMPI 4.1 (`/usr/bin/mpirun`). Test coverage:

| Suite | Tests |
|-------|-------|
| elastostatic (SPOOLES, N-LS) | 2 |
| elastostatic PCG / PCG.K | 2 |
| explicit dynamics + restart | 2 |
| split I/O elastostatic | 1 |
| 2D particle (PCG, MD, restart) | 3 |
| 2D particle periodic BC | 2 |
| 3D particle EAM / qsEAM | 2 |
| per-rank sub-tests (particle) | 8 |
| **Total** | **22** |

Run with:
```bash
./run_benchmarks.sh level.3
```

## Troubleshooting

### Wrong MPI picked up (conda MPICH)

Symptom: `ldd build/bin/tahoe` shows a conda path like
`libmpi.so => /home/user/miniforge3/lib/libmpi.so`.

Fix: wipe the CMake cache so stale `MPI_CXX_LINK_FLAGS` are discarded:
```bash
cmake -B build --fresh -DTAHOE_MPI=ON
```
The `--fresh` flag removes `build/CMakeCache.txt` entirely, forcing
re-detection of system MPI.

### `mpirun: command not found`

The system `mpirun` is at `/usr/bin/mpirun` (not in PATH if conda shadows it).
Either use the full path or install OpenMPI:
```bash
sudo apt-get install libopenmpi-dev   # installs /usr/bin/mpirun
```

### Link error: `curl@CURL_OPENSSL_4`

This occurs when conda's rpath (`-Wl,-rpath,/home/user/miniforge3/lib`) is
baked into `MPI_CXX_LINK_FLAGS` from a previous stale cache. The fix is the
same as above: `cmake -B build --fresh -DTAHOE_MPI=ON`.
