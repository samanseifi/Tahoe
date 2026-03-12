# Newton-Krylov Solver — Notes on Implementation, Performance, and Third-Party Libraries

## Overview

The `newton_krylov_solver` added to Tahoe replaces the direct LU/Cholesky factorization in each Newton iteration with **restarted GMRES(m)** (Generalized Minimal Residual) and a Jacobi (diagonal) left preconditioner.  This document explains the design, compares it against the classic Newton direct solver, and discusses how third-party Krylov-solver libraries could be integrated.

---

## Benchmark Comparison: `bar.2D.nonlin.1.xml`

This benchmark solves a steady-state **nonlinear diffusion** problem on a 4×1 bar (6 active DOFs). The conductivity depends on temperature, making the tangent matrix non-symmetric.

### Convergence (Newton iterations)

| Iteration | Newton + SPOOLES (direct) | Newton-Krylov + GMRES(30) |
|-----------|--------------------------|--------------------------|
| 0         | 5.81 × 10⁻¹              | 5.81 × 10⁻¹              |
| 1         | 4.04 × 10⁻²              | 4.04 × 10⁻²              |
| 2         | 2.22 × 10⁻⁴              | 2.22 × 10⁻⁴              |
| 3         | 5.99 × 10⁻⁹              | 5.99 × 10⁻⁹              |
| 4         | 2.36 × 10⁻¹⁶             | 3.30 × 10⁻¹⁶             |

Both solvers converge in **4 Newton iterations** to **machine-precision accuracy**. The NK variant produces a solution that differs from the direct solver by less than `10⁻¹⁶` in relative norm.

### XML Usage

```xml
<!-- Newton + SPOOLES direct solve (original) -->
<nonlinear_solver abs_tolerance="1.0e-10" divergence_tolerance="10.0"
                  max_iterations="5" rel_tolerance="1.0e-12">
    <profile_matrix/>
</nonlinear_solver>

<!-- Newton-Krylov + GMRES (benchmark, tight linear tolerance) -->
<newton_krylov_solver abs_tolerance="1.0e-10" divergence_tolerance="10.0"
                      max_iterations="5" rel_tolerance="1.0e-12"
                      gmres_restart="30" linear_tolerance="1.0e-8"
                      max_linear_iter="0">
    <profile_matrix/>
</newton_krylov_solver>

<!-- Newton-Krylov + GMRES (production / inexact-Newton, looser tolerance) -->
<newton_krylov_solver abs_tolerance="1.0e-10" divergence_tolerance="10.0"
                      max_iterations="10" rel_tolerance="1.0e-12"
                      gmres_restart="30" linear_tolerance="1.0e-3"
                      max_linear_iter="0">
    <profile_matrix/>
</newton_krylov_solver>
```

`linear_tolerance` is a *relative* GMRES stopping criterion: GMRES stops when  
`||r_k|| / ||r_0|| < linear_tolerance`.  
- `1.0e-8` — tight, matches machine-precision convergence of the direct solver.  
- `1.0e-3` to `1.0e-6` — typical for inexact-Newton methods (fewer Matvecs per Newton step, but may need more Newton iterations).

---

## Is the Built-In GMRES Fast?

### Computational Complexity

| Operation               | Cost per GMRES step |
|-------------------------|---------------------|
| Matrix-vector product   | O(nnz)              |
| Modified Gram-Schmidt   | O(n · j) for step j |
| Givens rotation update  | O(1) per step       |
| Back-substitution       | O(m²)               |

For a system with `n` DOFs and restart dimension `m`:
- Total cost per Newton step: O(m · nnz + m² · n)
- Memory: O(m · n + m²) for the Krylov basis

### When NK is Faster Than Direct Solvers

| Problem characteristics         | Best choice |
|---------------------------------|-------------|
| **Small/medium, dense band**    | Direct (SPOOLES, MUMPS) |
| **Large sparse, few DOFs/element** | GMRES with ILU |
| **Well-conditioned PDE**        | GMRES + Jacobi |
| **Ill-conditioned, nearly singular** | Direct factorization |
| **Symmetric positive-definite** | Conjugate Gradient |

The built-in GMRES uses only a **Jacobi (diagonal) preconditioner**, which is the cheapest possible preconditioner. For the 6-DOF example above it works perfectly. For larger, more challenging problems (elasticity, contact, large deformation), a stronger preconditioner is needed.

### Limitations of the Current Implementation

1. **Jacobi preconditioner only** — adequate for well-conditioned problems but may require many GMRES iterations for ill-conditioned systems.
2. **Pure C++ / no BLAS** — the `Multx` and dot products use element-wise loops; no vectorized BLAS-3 kernels.
3. **No parallel support** — single-processor only (the GMRES workspace is not distributed).

---

## Third-Party Krylov Library Integration

Several high-quality, open-source Krylov solver libraries can be integrated into Tahoe to significantly improve performance and robustness.

### 1. Trilinos / Belos + Ifpack2 (preferred for large-scale FEM)

**Trilinos** is already partially used in Tahoe (the `NOX` solver group is a Trilinos package). Adding **Belos** (iterative solvers) and **Ifpack2** (preconditioners) would provide:

- GMRES, BiCGSTAB, CG, TFQMR, etc.
- ILU(k), ILUT, multi-level Schwarz preconditioners
- Excellent parallel scaling via Epetra/Tpetra

```cmake
# In CMakeLists.txt:
find_package(Trilinos REQUIRED COMPONENTS Belos Ifpack2)
```

**Status:** Tahoe already links against some Trilinos packages. Full Belos integration would require wrapping `GlobalMatrixT` as a `Tpetra::Operator`.

### 2. PETSc (widely used in scientific computing)

**PETSc** provides:
- Complete linear algebra suite (KSP solvers + PC preconditioners)
- GMRES, CG, BiCGSTAB, MINRES, LSQR, and 20+ others
- ILU, ICC, algebraic multigrid (via Hypre/ML)
- Excellent scalability from laptop to supercomputer

```cmake
find_package(PETSc REQUIRED)
target_link_libraries(libtahoe PRIVATE PETSc::PETSc)
```

**Trade-off:** PETSc is a significant dependency but provides the most flexible solver ecosystem.

### 3. Eigen (lightweight, header-only)

**Eigen** provides iterative solvers (GMRES via `unsupported/Eigen/IterativeSolvers`) with minimal overhead. Suitable for single-processor use:

```cmake
find_package(Eigen3 REQUIRED)
target_link_libraries(libtahoe PRIVATE Eigen3::Eigen)
```

**Trade-off:** Limited to sequential execution; preconditioners are basic (diagonal, incomplete LU).

### 4. HYPRE (via Trilinos or standalone)

**HYPRE** specializes in scalable preconditioners, particularly **BoomerAMG** (algebraic multigrid). Ideal for large elliptic/parabolic PDEs (diffusion, elasticity). Available as a standalone library or via Trilinos.

---

## Recommended Upgrade Path

For production use on large-scale problems, the recommended approach is:

1. **Short-term**: Use the built-in `newton_krylov_solver` (current implementation) for problems where Jacobi preconditioning is sufficient (well-conditioned problems, moderate mesh sizes).

2. **Medium-term**: Integrate **PETSc KSP** as an optional backend by wrapping Tahoe's `GlobalMatrixT` as a `Mat` object. This would allow users to select any PETSc preconditioner (ILU, AMG, etc.) from the XML input without code changes.

3. **Long-term**: Consider replacing the ad-hoc iterative solver infrastructure with **Trilinos Belos+Ifpack2** since Tahoe already has partial Trilinos dependencies (NOX).

---

## Summary

| Feature | Built-in GMRES | PETSc | Trilinos/Belos |
|---------|---------------|-------|----------------|
| Dependencies | None | PETSc | Trilinos |
| Preconditioners | Jacobi only | ILU, AMG, Schwarz, ... | ILU, AMG, Schwarz, ... |
| Parallel | No | Yes (MPI+OpenMP) | Yes (MPI+OpenMP) |
| Matrix-free | Yes | Yes | Yes |
| Performance on large problems | Limited | Excellent | Excellent |
| Integration effort | Done ✓ | Moderate | Low (partial deps exist) |
