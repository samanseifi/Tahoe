# Tahoe Explicit Solver Performance Roadmap

## Current State (Issue11)

**Baseline:** Legacy `updated_lagrangian` with SPOOLES solver, single-threaded, no optimization flags.

**Achieved:** 9.3x speedup (20k elements, 10k steps, 6 OpenMP threads, -O3 -march=native)

| Configuration | Time (ms) | Speedup |
|--------------|-----------|---------|
| Legacy 4-IP | 173,271 | 1.0x |
| Legacy 1-IP | 51,348 | 3.4x |
| Vectorized 4-IP, 6 threads | 27,330 | 6.3x |
| Vectorized 1-IP, 6 threads | 18,534 | 9.3x |

## Completed Optimizations

1. **Diagonal mass auto-override** — bypasses sparse solver for lumped mass
2. **MVSIZ SoA batched element loop** — 128-element blocks for SIMD
3. **Single-pass Jacobian** — reference Jacobian + F push-forward, no double kernel call
4. **Virtual dispatch fix** — override RHSDriver (virtual) not ElementRHSDriver (non-virtual)
5. **OpenMP parallelism** — `#pragma omp parallel for` on batch loop
6. **Flat connectivity** — pre-computed array, no ElementCard virtual calls in gather

## Remaining Optimizations (CPU)

### Tier 1: Low effort, high impact

**A. Persistent force array** (1 hour)
- Allocate `dArray2DT force` once as a class member, not every time step
- Eliminates malloc + memset of (numnod × ndof) doubles per step
- Expected gain: 1.1-1.2x (significant for large meshes)

**B. Template specialization** (1 day)
- Template `BatchedInternalForce<NEN, NIP, NSD>()` on compile-time constants
- Compiler unrolls all inner loops, inlines shape derivatives as constants
- Eliminates virtual kernel dispatch per IP
- Expected gain: 1.5-2x

**C. Pre-computed reference shape derivatives** (2 hours)
- dN/dxi at each Gauss point are constants (same for all elements)
- Compute once at init, store as flat arrays
- Currently recomputed every time step inside the kernel
- Expected gain: 1.1-1.3x

**D. Profile-guided optimization (PGO) + link-time optimization (LTO)** (1 hour)
- `cmake -DCMAKE_CXX_FLAGS="-O3 -march=native -flto -fprofile-generate"`
- Run benchmark, then rebuild with `-fprofile-use`
- Expected gain: 1.2-1.5x (free, no code changes)

### Tier 2: Medium effort, medium impact

**E. Hourglass control for reduced integration** (2-3 days)
- Viscous hourglass stabilization for 1-point Q4 and Hex8
- Required for practical crash/impact simulations
- Not a speedup per se, but enables 1-IP elements which are 4-8x faster
- Formula: F_hg = rho * c * h * gamma_hg * (hourglass velocity modes)

**F. CFL-based adaptive time step** (1 day)
- dt = min over all elements of (h / c) where h = characteristic length, c = wave speed
- Currently using fixed dt — user must guess the stable time step
- Not a speedup, but prevents crashes and enables automatic time step selection
- Nodal formulation: dt_node = sqrt(M_node / K_node)

**G. Bulk viscosity for shock capturing** (1 day)
- Linear: Q1 = b1 * rho * c * h * |div(v)|
- Quadratic: Q2 = b2 * rho * h^2 * (div(v))^2
- Added to pressure in material stress computation
- Required for any problem with shocks or impacts

**H. Direct RHS scatter with pre-computed equation numbers** (half day)
- Pre-compute equation numbers per element at init (like flat connectivity)
- Scatter directly to global RHS pointer, bypassing AssembleRHS chain
- Expected gain: 1.1-1.2x (assembly is currently ~1.4% but grows with fewer IPs)

### Tier 3: High effort, high impact

**I. GPU offload (CUDA/HIP)** (2-3 weeks)
- Port inner MVSIZ loops to GPU kernels
- SoA layout is already GPU-friendly (coalesced memory access)
- Each GPU thread processes one element
- Gather on CPU → memcpy to GPU → compute → memcpy back → scatter on CPU
- For large meshes (100k+ elements): 50-200x over CPU
- Alternatively: unified memory to avoid explicit memcpy

**J. Full bypass of Tahoe infrastructure** (1 week)
- Bypass ContinuumElementT::RHSDriver (traction BCs handled separately)
- Bypass ElementSupport coordinate access (use raw pointers)
- Bypass dArray2DT operator() (use pointer arithmetic)
- Essentially: the explicit element becomes its own mini-solver
- Risk: loses compatibility with Tahoe features (output, restart, etc.)

**K. Multi-rate time stepping / subcycling** (1 week)
- Different dt for different element groups
- Small, stiff elements use smaller dt; large, soft elements use larger dt
- Reduces total step count by 2-5x for mixed meshes
- Complex: requires interpolation at interface between fast/slow groups

## Additional Material Support

### Batch materials to implement

| Material | Effort | Notes |
|----------|--------|-------|
| Ogden (hyperelastic) | 1 day | Eigenvalue decomposition per element, still vectorizable |
| J2 plasticity (radial return) | 2 days | History variables in SoA, branch divergence in SIMD (~20% penalty) |
| Mooney-Rivlin | half day | Similar to Neo-Hookean, different strain energy |
| Viscoelastic (Prony series) | 1 day | History: internal stress variables per Maxwell branch |
| Johnson-Cook (rate-dependent plasticity) | 2 days | Temperature + strain rate coupling |
| Equation of State (Mie-Gruneisen) | 1 day | Pressure-volume, needed for fluid/shock problems |
| Gurson (porous plasticity) | 3 days | Void fraction evolution, complex yield surface |

### Kernel topologies to implement

| Kernel | Effort | Notes |
|--------|--------|-------|
| Q1P0 (mean dilatation) | 2 days | Separate volumetric/deviatoric, extra mean pressure DOF |
| Hex8 1-point + hourglass | 2 days | Fastest 3D element for crash |
| Tet4 | half day | 1-IP linear tetrahedron, constant strain |
| Tet10 | 1 day | Quadratic tetrahedron, 4 IPs |
| Hex20/Hex27 | 2 days | Quadratic hexahedron, 27 IPs |
| Shell (MITC4) | 1 week | Thin structures, 5 DOF/node |
| Beam (Timoshenko) | 3 days | 1D structural elements |

## Contact and Multi-physics

### Contact (independent of explicit element optimization)

| Feature | Effort | Priority |
|---------|--------|----------|
| Bucket sort spatial search | 2 days | Critical for >1000 contact pairs |
| Surface-to-surface penalty | 3 days | More robust than node-to-segment |
| Coulomb friction | 1 day | Required for forming/crash |
| Tied interface | 1 day | Multi-material bonding |
| Self-contact | 2 days | Sheet metal forming |
| Mortar contact (for implicit) | 1 week | Higher accuracy, quasi-static |

### Multi-physics coupling

| Feature | Status | Notes |
|---------|--------|-------|
| Electromechanical (Maxwell stress) | Working | SimoQ1P0 + DiffusionElementT |
| Phase-field fracture | Working | PhaseFieldElementT |
| Surface tension | Working | SimoQ1P0_Surface |
| Thermal coupling | Partial | Heat generation exists, no thermal solve |
| Fluid (hydrostatic cavity) | Not started | Needed for HASEL actuators |
| FSI (ALE) | Not started | Fluid-structure interaction |

## Estimated Total Performance Envelope

| Configuration | Speedup vs baseline | Notes |
|--------------|-------------------|-------|
| Current (Issue11) | 9.3x | MVSIZ + OpenMP + flat gather |
| + Template specialization | 14-19x | Compile-time element type |
| + PGO/LTO | 17-28x | Free compiler optimization |
| + 1-IP hourglass | 70-110x | 4x fewer IPs, production-ready |
| + GPU | 500-1000x | Requires CUDA port |
| OpenRadioss (estimated) | 50-100x | Hand-tuned Fortran, decades of optimization |

## Build Configuration for Maximum Performance

```bash
# Release build with full optimization
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -DNDEBUG -flto" \
    -DTAHOE_OPENMP=ON

# Profile-guided optimization (two-pass)
# Pass 1: generate profile
cmake -DCMAKE_CXX_FLAGS="-O3 -march=native -fprofile-generate -flto" ..
make -j && ./tahoe -f benchmark.xml
# Pass 2: use profile
cmake -DCMAKE_CXX_FLAGS="-O3 -march=native -fprofile-use -flto" ..
make -j

# Run with optimal thread count (physical cores, not hyperthreads)
OMP_NUM_THREADS=6 ./tahoe -f input.xml
```
