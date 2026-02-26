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
| `TAHOE_F2C` | `ON` | Fortran-to-C converter (ABAQUS UMAT support) |
| `TAHOE_DEV` | `ON` | Research/development element module |
| `TAHOE_MPI` | `OFF` | MPI parallelization (requires system MPI) |
| `TAHOE_SEACAS` | `OFF` | ExodusII mesh I/O (requires ACCESS/SEACAS) |

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

### `spooles/` — Sparse Direct Solver
Bundled SPOOLES library (C, no external dependencies). Provides multi-frontal sparse LU and Cholesky factorization used as the default direct solver.

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

---

## Input Format

Tahoe reads XML input files validated against `tahoe.xsd`. Parameters are structured hierarchically using the `ParameterListT` system; sub-lists map directly to solver, material, and element configuration blocks.

---

## License

See the `LICENSE` file. Tahoe was developed at Sandia National Laboratories under US government funding and is distributed as open-source software.
