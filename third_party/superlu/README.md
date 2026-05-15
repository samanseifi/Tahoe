# SuperLU 3.0 — Bundled Serial Sparse Direct Solver

SuperLU 3.0 (July 2005) is a high-performance sparse direct solver that uses
column-permuted LU factorisation with partial pivoting. It is a drop-in
alternative to SPOOLES for serial (single-node, single-thread) jobs.

Activate with `-DTAHOE_SUPERLU=ON`. Requires `TAHOE_F2C=ON` (default).
No system BLAS installation is needed — the bundled `cblas/` directory
provides all required dense kernels as f2c-generated C code.

---

## Build

```bash
cmake -B build -DTAHOE_SUPERLU=ON      # TAHOE_F2C=ON is the default
cmake --build build -j$(nproc)
```

CMake guard: `TAHOE_SUPERLU=ON` will fail at configure time if
`TAHOE_F2C=OFF` is also specified.

---

## XML Usage

In any Tahoe XML input file, replace the solver sub-element of
`<nonlinear_solver>` with:

```xml
<nonlinear_solver abs_tolerance="1e-10" rel_tolerance="1e-8"
                  max_iterations="10" divergence_tolerance="1000">
    <SuperLU_matrix/>
</nonlinear_solver>
```

Optional parameters (both have defaults; omit to accept defaults):

```xml
<SuperLU_matrix print_stat="false" refinement="NOREFINE"/>
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `print_stat` | bool | `false` | Print SuperLU timing/fill statistics after each factorisation |
| `refinement` | enum | `NOREFINE` | Iterative refinement: `NOREFINE`, `SINGLE`, `DOUBLE`, `EXTRA` |

---

## Source Layout

```
superlu/
├── CMakeLists.txt       — builds tahoe_superlu static library
├── src/                 — SuperLU 3.0 C sources (127 files)
│   ├── dsp_defs.h       — double-precision SuperLU public API
│   ├── supermatrix.h    — SuperMatrix type definition
│   ├── util.h           — utility macros and types
│   ├── d*.c             — double-precision kernels (dgssvx, dgstrf, etc.)
│   ├── s*.c             — single-precision (built but not linked into Tahoe)
│   ├── sp_*.c           — shared sparse utilities (sp_ienv, sp_coletree, …)
│   ├── superlu_timer.c  — platform timer
│   ├── colamd.c/h       — COLAMD column-minimum-degree ordering
│   └── mmd.c            — Minimum degree ordering
└── cblas/               — Bundled double-precision BLAS (18 files, f2c-generated)
    ├── dgemm.c  dgemv.c  dtrsv.c  dtrsm.c
    ├── daxpy.c  dcopy.c  dscal.c  ddot.c
    └── …
```

### CMake targets

| Target | Type | Linked into |
|--------|------|-------------|
| `tahoe_superlu` | STATIC | `libtahoe` (when `TAHOE_SUPERLU=ON`) |

`tahoe_superlu` links `tahoe_f2c` privately to get `f2c.h` for the CBLAS
sources; the f2c include path is not exposed to consumers.

### Source filtering

The CMakeLists.txt builds only the precision variants needed by Tahoe:

| Pattern | Action | Rationale |
|---------|--------|-----------|
| `s*.c` (but not `sp_*.c`, `superlu_*.c`) | excluded | single-precision kernels |
| `c*.c` (but not `colamd.c`) | excluded | single-complex kernels |
| `z*.c` | excluded | double-complex kernels |
| `dreadhb.c` | excluded | Harwell-Boeing reader (not used by Tahoe) |
| everything else | included | shared utilities + double-precision kernels |

---

## Where it comes from

Sources are copied from the
[NguyenLabJHU/Tahoe-FEM](https://github.com/NguyenLabJHU/Tahoe-FEM) fork,
which bundles SuperLU 3.0 (released July 2005 by the SuperLU team at
UC Berkeley / Xerox PARC / Lawrence Berkeley National Lab) together with
a companion `SuperLU_CBLAS` library derived from SuperLU-DIST 2.0.

The original SuperLU 3.0 distribution is available at:
  https://portal.nersc.gov/project/sparse/superlu/

License: BSD-style (see copyright headers in `src/*.c`).

---

## Benchmark

The `benchmark_XML/level.1/material.solid/3D/material.120/` directory
contains a verification test that exercises SuperLU end-to-end:

| Input file | Material | Steps | Result |
|------------|----------|-------|--------|
| `wlc_superlu.xml` | Bischoff-Arruda WLC (finite anisotropy, 3D hex) | 290 | PASS |

Run from the `material.120/` directory:

```bash
cd benchmark_XML/level.1/material.solid/3D/material.120
/path/to/build/bin/tahoe -f wlc_superlu.xml
```

Or via the benchmark harness:

```bash
./run_benchmarks.sh level.1
```

Typical timing on a laptop (single-element problem):

| Solver | 290-step total |
|--------|---------------|
| SuperLU 3.0 | ~0.35 s |
| SPOOLES (serial) | ~0.40 s |

---

## Compiler flags

The CMakeLists.txt suppresses several C89/C90 warnings that are
harmless but noisy with modern GCC:

```cmake
-Wno-implicit-function-declaration
-Wno-implicit-int
-Wno-int-conversion
-Wno-unused-result
-Wno-return-type
-Wno-incompatible-pointer-types
```
