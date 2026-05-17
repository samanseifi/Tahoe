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
| `TAHOE_METIS` | `OFF` | System METIS 5 graph partitioner — improves MPI domain decomposition quality (minimises edge cuts / communication volume). Requires `TAHOE_MPI=ON`. Install: `sudo apt-get install libmetis-dev`. When enabled, METIS is used automatically during decomposition; disable at runtime with `-no_metis`. |
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

### `third_party/` — Bundled C libraries (built from source)

All six bundled C libraries live under `third_party/`.  None has external
dependencies; CMake builds them from source alongside the main project.

| Subdirectory | Role | Enable flag |
|---|---|---|
| [`third_party/spooles/`](third_party/spooles/README.md)       | SPOOLES 2.2 serial sparse direct solver (default) | `TAHOE_SPOOLES=ON` |
| [`third_party/spoolesMT/`](third_party/spoolesMT/README.md)   | POSIX-threads parallel LU on top of SPOOLES        | `TAHOE_SPOOLES_MT=ON` |
| [`third_party/spoolesMPI/`](third_party/spoolesMPI/README.md) | MPI distributed-memory LU on top of SPOOLES        | `TAHOE_MPI=ON` |
| [`third_party/superlu/`](third_party/superlu/README.md)       | SuperLU 3.0 serial sparse direct solver            | `TAHOE_SUPERLU=ON` |
| [`third_party/f2c/`](third_party/f2c/README.md)               | Fortran-to-C runtime (ABAQUS UMAT support)         | `TAHOE_F2C=ON` |
| [`third_party/expat/`](third_party/expat/README.md)           | XML parser for Tahoe's input format                | `TAHOE_EXPAT=ON` |

### MUMPS (system library, *not* bundled)
Wrapper around the system MUMPS direct solver (`libmumps-dev`). Two variants: serial (`MUMPSMatrixT` — uses `MPI_COMM_SELF`, no `mpirun` needed) and MPI distributed (`MUMPSMatrixT_mpi` — uses `MPI_COMM_WORLD`, requires `-DTAHOE_MPI=ON`). Enable with `-DTAHOE_MUMPS=ON`. Source: [`tahoe/src/primitives/globalmatrix/MUMPS/`](tahoe/src/primitives/globalmatrix/MUMPS/).

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

#### Benchmark Status (May 2026, default build — `cmake -B build`)

Use `run_benchmarks.sh` at the repo root to reproduce. Level 3 requires `-DTAHOE_MPI=ON` and system OpenMPI (`/usr/bin/mpirun`); the script auto-detects system `mpirun` and runs level.3 with 4 MPI ranks. The numbers below are the no-MPI default build except for level.3, which is reproduced from the MPI build.

| Level | PASS | FAIL/CRASH | SKIP | Notes |
|-------|------|------------|------|-------|
| level.0 | **172** | 14 | — | Core physics suite (up from 155 after issue #37 cleanup) |
| level.1 | **100** | 3 | — | Extended element tests (includes WLC + SuperLU benchmark) |
| level.2 | **47** | 0 | — | Additional verification (up from 39 — clean) |
| level.3 | **22** | 0 | — | MPI parallel (4 ranks): elastostatic, explicit dynamics, particle MD, PCG, periodic BC |
| **Total** | **341** | **17** | | up from 322/35 in Feb 2026 |

##### Recent improvements (issue #37, PR #38)
- Enabled `ADHESION_ELEMENT` in `tahoe/config/ElementsConfig.h` → 7 adhesion + Cauchy-Born tests recovered.
- Removed redundant `output_format="ExodusII"` from `diffusion/heat.0.xml` → format match restored.
- Reverted BC drift in `bar.2D.lin.xml` (`-8.0` → `-1.0`, regression from commit 48990950) → reference match restored.

##### Remaining failure categories

| Category | Count | Root cause |
|----------|-------|------------|
| **A** — Bridging element runtime crash | 9 | `BRIDGING_ELEMENT` flag compiles cleanly under the modernised build but crashes during atomistic-continuum coupling initialisation; flag kept off pending a separate investigation (sub-issue of #37). |
| **B** — Adhesion stress regression | 2 | `adhesion.{2,3}.xml` produce s11/s22 stresses that differ from reference by O(0.1); displacement field matches to ~1e-4. Numerical regression in adhesion element code path. |
| **C** — Output channel suppression | 2 | `inputoutput/square.xml` produces no `.io*.run` files; output system fails to register the configured I/O channels. |
| **D** — 5-node Q4 unsupported | 1 | `CSE.2.xml` has a bulk-element block with `nen=5` (Q4 + phantom-node CSE side); `QuadT::EvaluateShapeFunctions` only handles {1, 4, 8, 9} nodes. |
| **E** — Unregistered material | 3 | `mat.2.a.xml` and `mat.11.a.xml` (2D + 3D) reference materials not currently registered in the factory. |

Each remaining failure is tracked as a sub-investigation under issue #37; none are config-level fixes.

---

## Solvers

Tahoe ships five sparse direct solvers. The bundled solvers (SPOOLES, SuperLU) have no external library dependencies. The system solvers (MUMPS, SPOOLES-MPI) require system packages.

| Solver | CMake flag | Parallelism | XML element | When to use |
|--------|-----------|-------------|-------------|-------------|
| SPOOLES (default) | `TAHOE_SPOOLES=ON` | 1 thread | `<SPOOLES_matrix/>` | Default; development and small models |
| SuperLU 3.0 | `TAHOE_SUPERLU=ON` | 1 thread | `<SuperLU_matrix/>` | Serial alternative with partial pivoting; often faster than SPOOLES on medium models |
| MUMPS (serial) | `TAHOE_MUMPS=ON` | 1 thread | `<MUMPS_matrix/>` | System MUMPS on a single process; requires `libmumps-dev` |
| MUMPS (MPI) | `TAHOE_MUMPS=ON` + `TAHOE_MPI=ON` | N MPI ranks | `<MUMPS_MPI_matrix/>` | Distributed MUMPS across MPI ranks; same install, run with `mpirun -decomp_method -0` |
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

Verified: WLC finite-anisotropy benchmark (290 Newton steps, single hex element) completes in ~0.35 s with SuperLU vs ~0.40 s with SPOOLES. See [`third_party/superlu/README.md`](third_party/superlu/README.md).

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

<!-- MPI distributed — run with /usr/bin/mpirun and -decomp_method -0: -->
<MUMPS_MPI_matrix message_level="silent" always_symmetric="false"/>

<!-- message_level options: silent | errors | verbose -->
```

Both variants use AMD fill-reducing ordering (`icntl[6]=0`) and 200% workspace headroom (`icntl[13]=100`).

**How MUMPS-MPI works** (`icntl[17]=3`, distributed assembled input): each MPI rank calls `GenerateRCV()` to extract its local COO triplets with global indices, then supplies them to MUMPS directly. MUMPS handles the global assembly and distributed factorization internally. After the solve, the full solution is gathered on rank 0 and broadcast to all ranks, which extract their local portion.

> **BLAS note**: MUMPS links against system BLAS (`libblas.so.3`). The build uses `-Wl,--exclude-libs,ALL` to prevent SuperLU's bundled CBLAS routines from shadowing MUMPS's BLAS calls — without this, MUMPS segfaults on large problems inside the frontal factorization kernel.

> **mpirun note**: always use the system OpenMPI launcher (`/usr/bin/mpirun`), not conda's MPICH. With system OpenMPI, worker ranks redirect their output to per-rank `console<N>` log files, avoiding N-fold duplicated messages. Conda's `mpirun` bypasses rank detection and duplicates all output.

> **Domain decomposition**: `<MUMPS_MPI_matrix/>` requires the same domain-decomposition setup as SPOOLES-MPI — pass `-decomp_method -0` on the command line (or include it in a batch file). This generates per-rank geometry partitions automatically at startup.

Verified: WLC benchmark (290 Newton steps) and dielectric elastomer benchmark (18832 DOFs, mixed-physics) complete correctly with both `<MUMPS_matrix/>` and `<MUMPS_MPI_matrix/>`. See [`benchmark_XML/level.4/README.md`](benchmark_XML/level.4/README.md) for a full solver comparison.

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

### METIS — graph partitioner for MPI domain decomposition

METIS improves the quality of the mesh partition passed to MPI ranks. Without it, Tahoe uses a built-in recursive bisection which can be unbalanced on unstructured meshes. METIS minimises edge cuts (fewer shared boundary nodes → less inter-process communication per Newton iteration).

```bash
sudo apt-get install libmetis-dev

cmake -B build -DTAHOE_MPI=ON -DTAHOE_METIS=ON
cmake --build build -j$(nproc)

# METIS is used automatically — disable at runtime with -no_metis:
/usr/bin/mpirun -np 4 ./build/bin/tahoe -f input.xml -decomp_method -0
# force built-in bisection instead:
/usr/bin/mpirun -np 4 ./build/bin/tahoe -f input.xml -decomp_method -0 -no_metis
```

METIS partitioning is invoked during the decomposition phase (before the first time step). It calls `METIS_PartGraphKway` with vertex weights set to the nodal cost reported by each element group, producing load-balanced partitions with minimised cut edges. The runtime flag `-no_metis` falls back to the built-in bisection without rebuilding.

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
| February–March 2026 | Saman Seifi (Boston University) | **Build-system modernisation**: replaced Makefile/`tahoe-manager` with CMake (root, `toolbox/`, `tahoe/`, plus `expat/`, `spooles/`, `f2c/`, `superlu/`, `contrib/`, `benchmark_XML/CMakeLists.txt`); bundled libraries built from source (no pre-built `.a` files); GCC 11 compatibility fixes (C++11 two-phase lookup, `-fpermissive` for legacy templates, `#undef min/max/abs` to suppress STL macro pollution from `f2c.h`); `Environment.h` `__GCC_4__` mapping for GCC ≥ 5; rewrote broken `inc/` symlink trees against real `src/` directories. **CI/CD**: Google Test unit test suite (36 tests), GitHub Actions pipelines (serial / MUMPS / MPI), ExodusII/SEACAS via system packages (`libexodusii-dev`, `libnetcdf-dev`); fix missing `return` in `PotentialT::MeanEnergy`. **Solvers**: SPOOLES multithreaded (`TAHOE_SPOOLES_MT`, pthreads); SPOOLES distributed via bundled `spoolesMPI` (`TAHOE_MPI`) with 22/22 level.3 MPI benchmarks passing; SuperLU 3.0 serial (`TAHOE_SUPERLU`) with bundled CBLAS — no system BLAS; system MUMPS serial (`MUMPS_matrix`, `MPI_COMM_SELF`) and MPI (`MUMPS_MPI_matrix`, `icntl[17]=3` distributed-assembled); BLAS symbol conflict fix (`-Wl,--exclude-libs,ALL`) enabling MUMPS on 18k+ DOF problems (MUMPS-MPI np=2 fastest tested at 5.2 s vs 7.5 s serial / 13.4 s SPOOLES on 18k-DOF dielectric elastomer); METIS 5 graph partitioner (`TAHOE_METIS`) — rewrote `GraphBaseT::Partition_METIS` against the METIS 5 API (`METIS_PartGraphKway`, replaces removed METIS 4 functions). **Benchmark coverage**: 322/357 across level.0–level.3 (failure modes documented in the benchmark status table). **Materials**: enable `FINITE_ANISOTROPY` family (Bischoff-Arruda WLC). |
| March–April 2026 | Saman Seifi (Boston University) | **Classic-Tahoe (implicit / `<updated_lagrangian>`) modernisation**.  **SimoQ1P0** (issues #2, #3, #4): full refactor to make the element 2D/3D-consistent (div-div terms `/3.0 → /NumSD()` in 4 locations; B-bar exponent `1/3 → 1.0/NumSD()`) and to support coupled electro-mechanical-fracture physics — Maxwell stress and electrical tangent via `s_electric_ij` / `c_electrical_ijkl`, `TensorTransformT::PushForward` for Voigt push-forward; phase-field degradation g(d)=(1−d)²+k applied to stress *and* tangent; surface-tension Q1P0 element (3D and 2D), multi-block element groups with internal-interface surface tension (issue #10).  **Quadratic Newton convergence**: replaced the approximate analytical Q1P0 tangent with a numerical tangent (forward finite differences of `ComputeInternalForce`) — Newton goes from linear to quadratic on near-incompressible problems (#4).  **Phase-field fracture element** (AT2 model) added with a standalone verification benchmark (#4); **three-way electro-mechano-fracture coupling** between SimoQ1P0, the phase-field element, and a dielectric field.  **Linear-dielectric material** (constant ε, no electrostriction) registered for use as the bulk dielectric law (#3).  **Newton–Krylov solver**: `NewtonKrylovSolver` with GMRES(m) and `MSRMatrixT::Multx` for matrix-vector products — registered in the `FEManagerT` solver list (PR #7).  **Static dielectric-elastomer cleanup**: removed obsolete staggered formulations and the `DEDiffusion` element; resolved `FSDielectricElastomerQ1P0ElastocapillaryT` ODR collision; fixed `R` array sizing in `FormKd` to match `fRHS`; made the electric-potential field optional in SimoQ1P0; sideset-driven surface-tension input.  **Implicit Tet4 + ANP-Tet4** (`BonetTet`, issues #27, #28): plain Tet4 in `<updated_lagrangian>` and Bonet-Burton 1998 / ELFORM=13 F-bar averaged variant on the classic implicit path; known limitation tracked in issue #29 (BonetTet residual is computed from F̄ but the inherited tangent is from F → Newton stalls under near-incompressibility, fix is the same numerical-tangent approach used for SimoQ1P0). |
| April–May 2026 | Saman Seifi (Boston University) | **Modernised explicit solver track** (`ExplicitElementT`) with batch-vectorised internal-force kernels, OpenMP parallelism, single-pass Jacobian, flat connectivity arrays, mass scaling (fixed + adaptive), CFL-based time-step computation, viscous hourglass control, and material wave-speed propagation — issue #11.  Per-element Hex8/Q4 kernels (`Hex8Kernel`, `Q4Kernel`) with single-IP and full-IP variants; 29-116× speedup over the legacy explicit path on 3D and 2D benchmarks (level.5/explicit_benchmark).  **Batch finite-strain J2 plasticity** in the explicit element loop — issue #16.  **Tet4 element kernel** with Ji-transpose bug fix in all explicit kernels — issue #27.  **ANP-Tet4** (Bonet-Burton 1998 / LS-DYNA ELFORM=13 equivalent) F-bar averaged tetrahedron, both explicit and classic-Tahoe paths — issue #28.  **PenaltyContact3DT in explicit** — verified force-balance with central-difference, MVSIZ-batched striker loop — issue #19.  **Coulomb friction in PenaltyContact3DT** — regularised kinetic friction, retards sliding cube ~10 % at t=0.5 µs — issue #26.  **Contact-stack performance and stability fixes** — viscous damping (#31), explicit-element OpenMP auto-tune threshold (#32), parallel contact loop with thread-local workspaces (#33).  **Taylor bar benchmark** — issue #18; soft-metal demo, 60 µs validation, captures KE → plastic-work conversion.  **Hertz contact benchmark** — quarter-symmetry sphere indenter on elastic base, validates Hertz analytical to within 1-2 % on contact radius, peak pressure, and total load (level.5/hertz).  **Hertz Tet/Hex contact** (PR #35) — Tet4 indenter on Hex8 base (mixed-element contact via `Contact3DT::ConvertQuadToTri`), 1-3 % match to analytical with explicit ↔ implicit cross-validation < 0.5 %.  **Google Test suite expanded** to 50+ tests covering kernels, materials, contact-perf integration tests, and end-to-end benchmark verification. |
| May 2026 | Saman Seifi (Boston University) | **BonetTetT lagged-J̄ + FD tangent** (issue #29, PR #36) — `InitStep` caches J̄ at step start; `RHSDriver`/`LHSDriver` hold it frozen across all Newton iterations; `FormStiffness` builds the per-element 12×12 tangent by forward FD on the F̄ residual at frozen J̄. Both residual and tangent see the same `fJbarE` → consistent ∂R(u; J̄_step)/∂u and quadratic Newton. **15× iteration reduction** on `tet4_hyperelastic_anp.xml` (16 → 1 iter/step). Confirmed against Bonet & Burton 1998 (added to `ref/bonet1998.pdf`): F̄ formulation gives F̂_bar = F̂ and det(F̄) = J̄, equivalent to the iso/vol split in eqs (24)–(30) of the paper. Limitation documented in `BonetTetT.h`: cross-element ∂J̄/∂u coupling cannot fit in a per-element 12×12 block; severe κ ≫ μ remains future work (JFNK or analytical F-bar tangent à la Bonet/Marriott/Hassan 2001). **Failing-benchmark cleanup** (issue #37, PR #38) — enabled `ADHESION_ELEMENT` in `ElementsConfig.h`; dropped redundant `output_format="ExodusII"` from `diffusion/heat.0.xml`; reverted XML BC drift in `bar.2D.lin.xml`. Level.0 benchmark pass-rate: **172/186 PASS** (up from 155). Total across levels 0–3: **341 PASS / 17 FAIL** (up from 322/35). Remaining 17 failures are real Tahoe code regressions (bridging runtime crash, adhesion stress regression, output-channel suppression, 5-node Q4 unsupported, unregistered materials), tracked as sub-investigations under #37. |
