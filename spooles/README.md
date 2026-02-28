# SPOOLES 2.2 — Serial Sparse Direct Solver

Bundled copy of the SPOOLES 2.2 library from Rice University, built as a
static CMake target (`tahoe_spooles`). Provides multi-frontal sparse LU and
Cholesky factorisation for Tahoe's default single-process solver path.

## What SPOOLES provides

| Feature | Details |
|---------|---------|
| Factorisation | Multi-frontal sparse LU (asymmetric) and Cholesky (symmetric) |
| Ordering | Nested dissection, minimum degree, and custom permutations |
| Pivoting | Threshold pivoting for numerical stability |
| Interface | C API; wrapped in C++ by `SPOOLESMatrixT` |

## Integration with Tahoe

The CMake target is `tahoe_spooles`. Enabling it is controlled by:

```cmake
cmake -B build -DTAHOE_SPOOLES=ON   # ON by default
```

The C++ wrapper lives in the Tahoe source tree:

```
tahoe/src/primitives/globalmatrix/SPOOLES/
    SPOOLESMatrixT.h / .cpp    serial LU solver wrapper
```

`SolverT.cpp` selects `SPOOLESMatrixT` when `__SPOOLES__` is defined and
the problem runs on a single MPI rank (or MPI is not compiled in).

## Build targets

| CMake target | Library file | Purpose |
|---|---|---|
| `tahoe_spooles` | `libtahoe_spooles.a` | Serial SPOOLES |
| `tahoe_spooles_mt` | `libtahoe_spooles_mt.a` | Pthreads extension — see [`../spoolesMT/README.md`](../spoolesMT/README.md) |
| `tahoe_spooles_mpi` | `libtahoe_spooles_mpi.a` | MPI extension — see [`../spoolesMPI/README.md`](../spoolesMPI/README.md) |

## Source layout

```
spooles/src/
    A2/          Dense matrix operations (A2 object)
    BKL/         Block Krebs–Lehmer ordering
    BPG/         Bi-partite graph
    Chv/         Chevron frontal matrix
    ChvList/     Chevron list management
    ChvManager/  Memory manager for chevrons
    Coords/      Coordinate storage
    DSTree/      Domain/separator tree
    DenseMtx/    Dense matrix storage
    ETree/       Elimination tree
    FrontMtx/    Frontal matrix (key data structure)
    Graph/       Sparse graph connectivity
    ILUMtx/      Incomplete LU
    InpMtx/      Input matrix (triplet or CSC format)
    Lock/        Thread-locking primitives (serial: no-ops)
    MatMul/      Matrix multiplication utilities
    Misc/        Utility routines
    MSMD/        Multi-stage minimum degree ordering
    Network/     Network flow (for ordering)
    Perm/        Permutation vectors
    SymbFac/     Symbolic factorisation
    Tree/         General tree structure
    Utilities/   Sorting, hashing
    ZV/          Complex vector operations
    drivers/     High-level LU driver (library code — not test programs)
```

## Compile defines

| Define | Source | Effect |
|--------|--------|--------|
| `__SPOOLES__` | set by `tahoe_spooles` PUBLIC | Activates SPOOLES solver path in `SolverT.cpp` and `toolbox` |
| `__POSIX_THREADS__` | set when `TAHOE_SPOOLES_MT=ON` | Selects `pthread_mutex_t` layout in `Lock.h`; safe to add even for serial builds because serial drivers pass `NO_LOCK` |

## Upstream

SPOOLES 2.2, Clay Ashcraft & Cleve Ashcraft, Rice University.
Original distribution: http://www.caam.rice.edu/software/SPOOLES/
