# Unit Tests

Google Test suite for core Tahoe data structures, shape functions, and
material models. Tests are built automatically when `TAHOE_TESTS=ON`
(the default) and run with `ctest`.

## Running

```bash
# Build (tests are built automatically)
cmake -B build && cmake --build build -j$(nproc)

# Run all 36 tests
ctest --test-dir build

# Verbose output on failure
ctest --test-dir build --output-on-failure

# Run a single test suite
ctest --test-dir build -R "^NeoHookean"
```

## Test inventory (36 tests)

### `toolbox/test_dArrayT.cpp` — `dArrayT` (dynamic array)

| Test | Checks |
|------|--------|
| `DefaultConstructorIsEmpty` | Zero length after default construction |
| `ExplicitDimension` | Correct size after `Dimension(N)` |
| `ScalarAssignment` | All elements set by scalar `operator=` |
| `Magnitude` | L2 norm of a known vector |
| `OperatorPlusEquals` | In-place addition |
| `ScalarMultiply` | Scalar scale operator |
| `CopyAssignment` | Deep copy semantics |

### `toolbox/test_dMatrixT.cpp` — `dMatrixT` (dense matrix)

| Test | Checks |
|------|--------|
| `Det2x2` | Determinant of a 2×2 matrix |
| `Det3x3Identity` | det(I₃) = 1 |
| `Trace` | Sum of diagonal elements |
| `Inverse2x2` | A · A⁻¹ = I |
| `Symmetrize` | Symmetrisation operation |
| `ScalarAssignment` | All elements set to scalar |

### `toolbox/test_dSymMatrixT.cpp` — `dSymMatrixT` (symmetric matrix, Voigt storage)

| Test | Checks |
|------|--------|
| `Dimension3DStorageSize` | 6 Voigt components for 3D |
| `Dimension2DStorageSize` | 3 Voigt components for 2D |
| `Trace3D` | Trace of 3D symmetric matrix |
| `Det3DDiagonal` | Determinant of diagonal 3D matrix |
| `DeviatoricTraceZero` | tr(dev(A)) = 0 |
| `Eigenvalues3DDiagonal` | Eigenvalues of diagonal matrix = diagonal entries |

### `geometry/test_QuadT.cpp` — `QuadT` (bilinear quadrilateral shape functions)

| Test | Checks |
|------|--------|
| `PartitionOfUnityAtCenter` | Σ Nᵢ(0,0) = 1 |
| `PartitionOfUnityAtGaussPoint` | Σ Nᵢ at Gauss points = 1 |
| `NodalRecoveryNode0` | N₀ = 1 at node 0, 0 at all others |
| `NodalRecoveryNode1` | N₁ = 1 at node 1, 0 at all others |
| `ShapeFunctionsNonNegativeInterior` | Nᵢ ≥ 0 inside the element |

### `materials/test_IsotropicT.cpp` — `IsotropicT` (isotropic elastic moduli)

| Test | Checks |
|------|--------|
| `DefaultConstructorYieldsZeroModuli` | μ = κ = 0 after default construction |
| `ShearModulusFromENu` | μ = E / 2(1+ν) |
| `LameModulusFromENu` | λ = Eν / (1+ν)(1-2ν) |
| `YoungFromMuKappa` | E recovered from (μ, κ) |
| `NegativeYoungThrows` | Exception for E < 0 |
| `PoissonOutOfRangeThrows` | Exception for ν ∉ (-1, 0.5) |

### `materials/test_NeoHookean.cpp` — Neo-Hookean hyperelastic material

| Test | Checks |
|------|--------|
| `EnergyAtIdentity` | W(I) = 0 (stress-free reference state) |
| `DevStressZeroAtIdentity` | τ_dev(I) = 0 |
| `DevStressIsDeviatoric` | tr(τ_dev) = 0 for arbitrary deformation |
| `MeanStressAtJ2` | Hydrostatic pressure at J=2 matches analytical κ·ln(2) |
| `MeanStressZeroAtJ1` | Zero pressure at J=1 |
| `DevModDiagonalPositive` | Positive diagonal of deviatoric tangent modulus |

## Adding tests

Each subdirectory maps to one CMake target in `tests/CMakeLists.txt`. To add
a new test:

1. Create `tests/<category>/test_MyClass.cpp` with `#include <gtest/gtest.h>`
2. Register it in `tests/CMakeLists.txt` using the `tahoe_gtest` helper macro:

```cmake
tahoe_gtest(test_MyClass toolbox/test_MyClass.cpp)
```

For tests that need `libtahoe` (not just toolbox), use `tahoe_libtahoe_gtest`
if it has been defined, or add `libtahoe` to `target_link_libraries` manually.
