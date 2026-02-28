# SPOOLES-MT — Multithreaded Sparse Direct Solver

POSIX-thread parallel extension to SPOOLES 2.2. Parallelises the LU
factorisation and solve across multiple cores on a **single node** without
requiring MPI. Sourced from the NguyenLabJHU/Tahoe-FEM fork, which bundles
the original `spoolesMT` directory from the SPOOLES 2.2 distribution.

## Requirements

| Dependency | Notes |
|------------|-------|
| `TAHOE_SPOOLES=ON` | Serial SPOOLES must be built first |
| pthreads | Standard on Linux/macOS; CMake finds via `find_package(Threads)` |
| No MPI | Single-node only |

## Enabling

```bash
cmake -B build -DTAHOE_SPOOLES=ON -DTAHOE_SPOOLES_MT=ON
cmake --build build -j$(nproc)
```

## How it works

SPOOLES-MT adds a parallel factorisation layer (`MT.factorMT.c`,
`MT.solveMT.c`) on top of the serial SPOOLES frontal matrix data structures.
The driver (`LU_MT_driver`) accepts a thread count at initialisation and
distributes front-panel tasks across a pthread pool.

### Thread-safety note on `Lock.h`

SPOOLES's `Lock` object uses compile-time selection between a no-op stub and a
real `pthread_mutex_t`. The selector is `__POSIX_THREADS__`. When
`TAHOE_SPOOLES_MT=ON`, this define is propagated to **both** `tahoe_spooles`
(serial library) and `tahoe_spooles_mt`. This is safe: the serial LU driver
always passes `NO_LOCK (lockflag=0)` to `Lock_init`, so no mutex is ever
allocated in the serial code path.

## Source layout

```
spoolesMT/src/
    MT/
        MT.factorMT.c       parallel LU factorisation kernel
        MT.solveMT.c        parallel triangular solve
        MT.mvm.c            parallel matrix-vector multiply
        MT.QRfactorMT.c     parallel QR factorisation
        MT.QRsolveMT.c      parallel QR solve
        MT.spoolesMT.h      internal header
    drivers/
        LU_MT_driver.c      high-level entry point (init/factor/solve/free)
        LU_MT_driver.h      public API
        LU_MT_driver_int.h  internal driver header
        LU_MT_driver_init.c
        LU_MT_driver_factorize.c
        LU_MT_driver_solve.c
        LU_MT_driver_free.c
    SPOOLESMT.h             top-level include
```

## Tahoe wrapper

```
tahoe/src/primitives/globalmatrix/SPOOLES/
    SPOOLESMatrixT_MT.h / .cpp    pthreads LU wrapper
```

`SolverT.cpp` selects `SPOOLESMatrixT_MT` when both `__SPOOLES_MT__` and
`__POSIX_THREADS__` are defined.

## Compile defines

| Define | Set by | Effect |
|--------|--------|--------|
| `__SPOOLES_MT__` | `tahoe/CMakeLists.txt` (PUBLIC on libtahoe) | Activates `SPOOLESMatrixT_MT` in `SolverT.cpp` |
| `__POSIX_THREADS__` | `tahoe_spooles_mt` (PUBLIC) | Selects pthread layout in `Lock.h`; propagated to `tahoe_spooles` too |

## XML usage

In the solver block of your input XML, replace the default SPOOLES entry with:

```xml
<SPOOLES_MT_matrix num_threads="8" .../>
```

Set `num_threads` to the number of physical cores available. Over-subscription
(more threads than cores) typically reduces performance.

## Benchmark results

All 36 Google Test unit tests pass with `TAHOE_SPOOLES_MT=ON`. Level.0–2
benchmark results are identical to the serial build (no regressions).
