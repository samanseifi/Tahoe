<p align="center">
  <img src="tahoe_logo_v21.png" alt="Tahoe Logo" width="480"/>
</p>

**Tahoe** is an open-source, modular C++ finite element framework for solid mechanics, multiscale analysis, and advanced material modeling. Originally developed at Sandia National Laboratories, it provides a flexible research platform covering a broad range of element formulations, time integrators, and material models.

---

## Building

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
# Executables → build/bin/  |  Libraries → build/lib/
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `TAHOE_EXPAT` | `ON` | Bundled expat XML parser |
| `TAHOE_SPOOLES` | `ON` | Bundled SPOOLES sparse direct solver |
| `TAHOE_SPOOLES_MT` | `OFF` | SPOOLES multithreaded solver — POSIX-thread parallel LU factorisation on a single node; no MPI required. Requires `TAHOE_SPOOLES=ON`. Use `<SPOOLES_MT_matrix num_threads="N" .../>` in XML (N ≥ 2). |
| `TAHOE_SUPERLU` | `OFF` | Bundled SuperLU 3.0 serial sparse direct solver — high-performance LU with partial pivoting and optional iterative refinement. Requires `TAHOE_F2C=ON`. No system BLAS needed. Use `<SuperLU_matrix/>` in XML. |
| `TAHOE_MUMPS` | `OFF` | System MUMPS direct solver — links against system `libmumps-dev`. Two variants: `<MUMPS_matrix/>` (serial, uses `MPI_COMM_SELF`) and `<MUMPS_MPI_matrix/>` (distributed, requires `TAHOE_MPI=ON`). Install: `sudo apt-get install libmumps-dev`. |
| `TAHOE_F2C` | `ON` | Fortran-to-C converter (ABAQUS UMAT support) |
| `TAHOE_DEV` | `ON` | Research/development element module |
| `TAHOE_MPI` | `OFF` | MPI parallelization — requires system OpenMPI (`libopenmpi-dev`). Automatically builds the bundled `spoolesMPI` distributed solver. CMake prefers system wrappers (`/usr/bin/mpicxx`) over conda-installed MPI; override with `-DMPI_CXX_COMPILER=...`. Run with `mpirun -np N tahoe -f input.xml`. |
| `TAHOE_SEACAS` | `OFF` | ExodusII mesh I/O — auto-detects system packages or `ACCESS` tree (see below) |
| `TAHOE_TESTS` | `ON` | Build Google Test unit test suite |

#### Enabling ExodusII (SEACAS)

Two discovery modes are tried automatically when `-DTAHOE_SEACAS=ON`:

1. **System packages** (Ubuntu/Debian): install `libexodusii-dev libnetcdf-dev`, then:
   ```bash
   cmake -B build -DTAHOE_SEACAS=ON
   ```
2. **SEACAS/ACCESS tree**: set the `ACCESS` environment variable or pass `-DACCESS_PATH=...`.

---

## Repository Structure

### `tahoe/` — Core FEM Engine
The main analysis library and executable. Implements the primary element library, time integration, solvers, and I/O.

**Elements**
| Subdirectory | Description |
|---|---|
| `continuum/solid` | Small- and large-strain solid continuum elements (Q1, Q1P0, Simo enhanced, etc.) |
| `continuum/diffusion` | Diffusion and heat transfer elements |
| `cohesive_surface` | Cohesive zone models for fracture (iso- and anisotropic) |
| `contact` | Penalty and constraint-based contact elements |
| `bridging_scale` | Atomistic-to-continuum bridging scale elements |
| `adhesion` | Adhesive contact with traction-separation laws |
| `particle` | Discrete particle elements |
| `spring` | Truss and spring elements |
| `shape_functions` | Standard FEM and meshfree (EFG, reproducing kernel) shape functions |
| `constant_volume` | Constant-volume constraint elements |

**Integrators** — Static, explicit central difference, Verlet, trapezoidal, HHT-α, Newmark, Gear 6th-order, mixed

**Solvers** — Full Newton–Raphson, line search, DR (dynamic relaxation), preconditioned CG

---

### `development/` — Research Elements
Advanced and experimental formulations, enabled with `-DTAHOE_DEV=ON`. Several sub-modules are guarded by feature flags and compiled selectively.

| Module | Description |
|---|---|
| `APS_grad` | Gradient-enhanced plasticity (Aifantis-type) |
| `APS_grad_vector` | Vector-valued APS gradient formulation |
| `AG_element` | Alternative geometry element (quadratic, Andelfinger–Ramm) |
| `cohesive_surface` | Extended cohesive models (inelastic, ductile) |
| `contact` | Development contact elements |
| `DE` / `DEQ1P0` / `DEQ1P02D` | Dielectric elastomer elements (2D/3D, mixed pressure) |
| `DEQ1P0_Elastocapillary` | Coupled elastocapillary dielectric formulation |
| `Visco_DE` variants | Viscoelastic dielectric elastomer elements |
| `fluid_element` | Incompressible Navier–Stokes finite elements |
| `HuWashizu` | Hu–Washizu three-field mixed elements |
| `meshfree_grad_plast` | Meshfree gradient plasticity |
| `micromorphic` / `micromorphic2` | Higher-order micromorphic continuum (2D and 3D) |
| `micromorphic_curr_config` | Micromorphic formulation in current configuration |
| `multiscale` / `FEA` | FE² computational homogenization (multiscale) |
| `piezoelectric` | Piezoelectric coupled elements |
| `solid` | Additional solid element formulations |
| `solid_fluid_mix` | Solid–fluid mixture (poromechanics, Darcy) |
| `surface_CB` / `surface_CB_EAM` / `surface_CB_Si` | Surface Cauchy–Born atomistic-informed elements |
| `xfem` | Extended FEM for crack propagation |
| `optimization` | Sensitivity and topology optimization support |

---

### `toolbox/` — Utility Library
Low-level infrastructure shared by all modules.

| Subdirectory | Description |
|---|---|
| `abc/` | Core array, matrix, and linear algebra templates (`dArrayT`, `dMatrixT`, `nMatrixT`, etc.) |
| `dataio/` | Mesh I/O (ExodusII, HDF5/XDMF, text, VTK) and output management |
| `param_tree/` | XML-driven parameter input system (`ParameterListT`, `ParameterInterfaceT`) |
| `geometry/` | Geometric algorithms, bounding boxes, convex hull |
| `graph/` | Sparse graph connectivity and reordering |
| `search/` | Spatial search structures (grid manager, neighbor search) |
| `parallel/` | MPI communication wrappers |
| `C1functions/` | C¹-continuous interpolation functions (splines, piecewise linear) |
| `linkedlist/` | Linked list, binary tree, and map containers |
| `neighbors/` | Neighbor-list construction for particles/meshfree methods |
| `misc/` | String utilities, exception handling, timing |

---

### `spooles/` — Sparse Direct Solver (serial)
Bundled SPOOLES 2.2 library (C, no external dependencies). Provides multi-frontal sparse LU and Cholesky factorization; the default solver for serial and shared-memory builds. See [`spooles/README.md`](spooles/README.md).

### `spoolesMT/` — SPOOLES Multithreaded
POSIX-threads extension to SPOOLES. Parallelises the LU factorisation across cores on a single node without MPI. Enable with `-DTAHOE_SPOOLES_MT=ON`. See [`spoolesMT/README.md`](spoolesMT/README.md).

### `spoolesMPI/` — SPOOLES Distributed (MPI)
MPI extension to SPOOLES for distributed-memory parallel factorisation across multiple nodes. Enable with `-DTAHOE_MPI=ON`. See [`spoolesMPI/README.md`](spoolesMPI/README.md).

### `superlu/` — SuperLU 3.0 Sparse Direct Solver
Bundled SuperLU 3.0 library (C, no system BLAS dependency). Provides sparse LU factorisation with partial pivoting and optional iterative refinement; a fast drop-in alternative to SPOOLES for serial single-node jobs. Enable with `-DTAHOE_SUPERLU=ON`. See [`superlu/README.md`](superlu/README.md).

### MUMPS (system library)
Wrapper around the system MUMPS direct solver (`libmumps-dev`). Two variants: serial (`MUMPSMatrixT` — uses `MPI_COMM_SELF`, no `mpirun` needed) and MPI distributed (`MUMPSMatrixT_mpi` — uses `MPI_COMM_WORLD`, requires `-DTAHOE_MPI=ON`). Enable with `-DTAHOE_MUMPS=ON`. Source: [`tahoe/src/primitives/globalmatrix/MUMPS/`](tahoe/src/primitives/globalmatrix/MUMPS/).

### `f2c/` — Fortran-to-C Runtime
Enables ABAQUS UMAT material subroutines (originally written in Fortran) to be compiled and called from C++.

### `expat/` — XML Parser
Bundled expat library. Parses Tahoe's XML input format (validated against `tahoe.xsd`).

### `contrib/` — Pre/Post-Processing Tools
| Tool | Description |
|---|---|
| `MakeCSE` | Comparative Suite Environment — batch test runner |
| `translate` | Mesh format conversion (Patran, Abaqus → Tahoe XML) |
| `vtk` | VTK/ParaView output adaptor |
| `Dakota` | Integration with Dakota for optimization and UQ |
| `wrap` | Code-generation wrappers |
| `cubit` | Cubit mesh journal file support |

### `benchmark_XML/` — Regression Tests
18+ benchmark problems covering elastostatics, elastodynamics, diffusion, contact, cohesive fracture, particle methods, meshfree analysis, and time integrator verification. Run through `MakeCSE`.

See `benchmark_XML/ReadMe` for detailed run instructions. Quick start (from repo root):
```bash
# Run all levels and print a pass/fail summary
./run_benchmarks.sh

# Or run a single level
./run_benchmarks.sh level.0

# Manual run for one directory (from build root)
cd benchmark_XML/level.0
printf "run.batch\nquit\n" | ../../build/bin/tahoe   # run simulations
printf "run.batch\nquit\n" | ../../build/bin/compare  # compare vs reference
```

#### Benchmark Status (February 2026, MPI build with SEACAS enabled)

Use `run_benchmarks.sh` at the repo root to reproduce. Level 3 requires `-DTAHOE_MPI=ON` and system OpenMPI (`/usr/bin/mpirun`). The script auto-detects system `mpirun` and runs level.3 with 4 MPI ranks.

| Level | PASS | FAIL/CRASH | SKIP | Notes |
|-------|------|------------|------|-------|
| level.0 | **155** | 30 | — | Core physics suite |
| level.1 | **106** | 3 | — | Extended element tests (includes WLC + SuperLU benchmark) |
| level.2 | **39** | 2 | — | Additional verification |
| level.3 | **22** | 0 | — | MPI parallel (4 ranks): elastostatic, explicit dynamics, particle MD, PCG, periodic BC |
| **Total** | **322** | **35** | | |

##### Failure categories

| Category | Count | Root cause |
|----------|-------|------------|
| **A** — ExodusII format mismatch | ~16 | Tests produce `.exo` output; reference files use TahoeII `.run/.geo` format; `compare` cannot read ExodusII. Physics results are correct — regenerating references with SEACAS enabled would clear these. |
| **B** — Missing compiled features | ~11 | Bridging-scale element (`BRIDGING_ELEMENT` flag), surface Cauchy–Born (`surface_CB`), and CSE hex-shape tests require optional compile-time modules not built in the standard configuration. |
| **C** — SimoQ1P0 Voltage field | 5 | `SimoQ1P0` unconditionally searches for an electrical-field DOF at construction; pure-mechanical Q1P0 tests fail with *"Voltage field not found"*. |
| **D** — NaN in diffusion solver | 4 | `heat.1`, `heat.2`, `tsurf`, `heat.hyper.1` produce `-nan` temperatures; numerical issue in the diffusion element path. |
| **E** — Unregistered material | 2 | `small_strain_StVenant_DP_2D` is not registered in the material factory (`mat.11` tests). |

---

## Solvers

Tahoe ships five sparse direct solvers. The bundled solvers (SPOOLES, SuperLU) have no external library dependencies. The system solvers (MUMPS, SPOOLES-MPI) require system packages.

| Solver | CMake flag | Parallelism | XML element | When to use |
|--------|-----------|-------------|-------------|-------------|
| SPOOLES (default) | `TAHOE_SPOOLES=ON` | 1 thread | `<SPOOLES_matrix/>` | Default; development and small models |
| SuperLU 3.0 | `TAHOE_SUPERLU=ON` | 1 thread | `<SuperLU_matrix/>` | Serial alternative with partial pivoting; often faster than SPOOLES on medium models |
| MUMPS (serial) | `TAHOE_MUMPS=ON` | 1 thread | `<MUMPS_matrix/>` | System MUMPS on a single process; requires `libmumps-dev` |
| MUMPS (MPI) | `TAHOE_MUMPS=ON` + `TAHOE_MPI=ON` | N MPI ranks | `<MUMPS_MPI_matrix/>` | Distributed MUMPS across MPI ranks; same install, run with `mpirun` |
| SPOOLES-MT | `TAHOE_SPOOLES_MT=ON` | N pthreads | `<SPOOLES_MT_matrix num_threads="N"/>` | Multicore workstations; no MPI required |
| SPOOLES-MPI | `TAHOE_MPI=ON` | N MPI ranks | batch mode | HPC clusters or multi-node jobs |

### Serial SPOOLES (default)

No extra flags needed. Built automatically when `TAHOE_SPOOLES=ON`.

```bash
cmake -B build
cmake --build build -j$(nproc)
./build/bin/tahoe -f input.xml
```

### SuperLU 3.0 — serial, high-performance direct solver

Requires `TAHOE_F2C=ON` (on by default). No system BLAS needed — uses bundled CBLAS.

```bash
cmake -B build -DTAHOE_SUPERLU=ON
cmake --build build -j$(nproc)
./build/bin/tahoe -f input.xml
```

In the XML input, replace the solver block with:
```xml
<SuperLU_matrix/>
<!-- optional parameters: -->
<SuperLU_matrix print_stat="false" refinement="NOREFINE"/>
<!-- refinement options: NOREFINE | SINGLE | DOUBLE | EXTRA -->
```

Verified: WLC finite-anisotropy benchmark (290 Newton steps, single hex element) completes in ~0.35 s with SuperLU vs ~0.40 s with SPOOLES. See [`superlu/README.md`](superlu/README.md).

### MUMPS — system sparse direct solver (serial and MPI)

Requires `libmumps-dev`. Two variants share the same install and build flags:

- **Serial** (`<MUMPS_matrix/>`) — uses `MPI_COMM_SELF`; runs in-process without `mpirun`
- **MPI** (`<MUMPS_MPI_matrix/>`) — distributes factorization across ranks; requires `TAHOE_MPI=ON`

```bash
sudo apt-get install libmumps-dev

# Serial MUMPS only (recommended):
cmake -B build -DTAHOE_MUMPS=ON
cmake --build build -j$(nproc)
./build/bin/tahoe -f input.xml

# MUMPS + MPI (enables both variants; use system OpenMPI launcher):
cmake -B build -DTAHOE_MUMPS=ON -DTAHOE_MPI=ON
cmake --build build -j$(nproc)
/usr/bin/mpirun -np 4 ./build/bin/tahoe -f input.xml
```

In the XML input, replace the solver block with:
```xml
<!-- Serial (no mpirun needed) — recommended: -->
<MUMPS_matrix message_level="silent" always_symmetric="false"/>

<!-- MPI distributed (run with /usr/bin/mpirun, experimental): -->
<MUMPS_MPI_matrix message_level="silent" always_symmetric="false"/>

<!-- message_level options: silent | errors | verbose -->
```

Uses AMD fill-reducing ordering (`icntl[6]=0`) and 200% workspace headroom (`icntl[13]=100`).

> **BLAS note**: MUMPS links against system BLAS (`libblas.so.3`). The build uses `-Wl,--exclude-libs,ALL` to prevent SuperLU's bundled CBLAS routines from shadowing MUMPS's BLAS calls — without this, MUMPS segfaults on large problems inside the frontal factorization kernel.

> **mpirun note**: always use the system OpenMPI launcher (`/usr/bin/mpirun`), not conda's MPICH. With system OpenMPI, worker ranks redirect their output to per-rank `console<N>` log files, avoiding N-fold duplicated messages. Conda's `mpirun` bypasses rank detection and duplicates all output.

> **MUMPS_MPI_matrix note**: experimental — requires proper domain decomposition (batch mode, same requirements as SPOOLES-MPI). For large single-node problems, `<MUMPS_matrix/>` (serial MUMPS without mpirun) is the recommended choice.

Verified: WLC benchmark (290 Newton steps) and dielectric elastomer benchmark (18832 DOFs, mixed-physics) complete correctly with `<MUMPS_matrix/>`.

### SPOOLES-MT — shared-memory multithreaded

Requires pthreads (standard on Linux/macOS). No MPI installation needed.

```bash
cmake -B build -DTAHOE_SPOOLES_MT=ON
cmake --build build -j$(nproc)
./build/bin/tahoe -f input.xml   # thread count set inside the XML
```

In the XML input, replace the solver block with:
```xml
<SPOOLES_MT_matrix num_threads="8" .../>
```

### SPOOLES-MPI — distributed-memory MPI

Requires system OpenMPI (`libopenmpi-dev` on Debian/Ubuntu). CMake auto-detects `/usr/bin/mpicxx`; pass `-DMPI_CXX_COMPILER=/path/to/mpicxx` to override.

```bash
# Install system OpenMPI first (not conda — see note below)
sudo apt-get install libopenmpi-dev

cmake -B build -DTAHOE_MPI=ON
cmake --build build -j$(nproc)
mpirun -np 4 ./build/bin/tahoe -f input.xml
```

> **Conda note**: if your shell has conda activated, `mpirun` may point to conda's MPICH, which is ABI-incompatible with the system OpenMPI that Tahoe links against. The CMake scripts explicitly prefer `/usr/bin/mpicxx` over conda wrappers. After a fresh build, verify with `ldd build/bin/tahoe | grep mpi` — it should show `/lib/x86_64-linux-gnu/libmpi.so`, not a conda path. If you see a conda path, run `cmake -B build --fresh -DTAHOE_MPI=ON` to wipe the stale cache.

**Batch-mode invocation** (used by `run_benchmarks.sh` for level.3):
```bash
cd benchmark_XML/level.3/parallel
mpirun -np 4 /path/to/tahoe -f run.batch
```

The batch file format uses `@` as the first character (batch-mode marker), followed by option lines (`-run`, `-join_io`, `-decomp_method -0`, etc.) and XML filenames. Options accumulate as global state before each XML file is processed.

---

## Input Format

Tahoe reads XML input files validated against `tahoe.xsd`. Parameters are structured hierarchically using the `ParameterListT` system; sub-lists map directly to solver, material, and element configuration blocks.

---

## License

See the `LICENSE` file. Tahoe was developed at Sandia National Laboratories under US government funding and is distributed as open-source software.

---

## Revision History

| Date | Author | Notes |
|------|--------|-------|
| 2014 | Regents of the University of Colorado | Tahoe 2.1 release |
| February 2026 | Saman Seifi (Boston University) | CMake modernization; C++11 two-phase lookup fixes; compiler warning cleanup (`-fpermissive`); ExodusII/SEACAS enabled via system packages; Google Test unit test suite (36 tests); GitHub Actions CI/CD pipeline; fix missing `return` in `PotentialT::MeanEnergy`; SPOOLES multithreaded (pthreads) solver (`TAHOE_SPOOLES_MT`); MPI distributed solver via bundled spoolesMPI (`TAHOE_MPI`); 22/22 level.3 MPI benchmarks pass; SuperLU 3.0 serial solver (`TAHOE_SUPERLU`) with bundled CBLAS — no system BLAS required; enable `FINITE_ANISOTROPY` materials (Bischoff-Arruda WLC); system MUMPS solver — serial (`MUMPS_matrix`) and MPI distributed (`MUMPS_MPI_matrix`); fix BLAS symbol conflict (`-Wl,--exclude-libs,ALL`) enabling MUMPS on large problems (18k+ DOFs); 322/357 benchmarks pass |
