# ---------------------------------------------------------------------------
# SuperLU 3.0 — serial sparse direct solver (built from bundled source)
# BLAS routines provided by the bundled SuperLU_CBLAS (f2c-generated C code
# from SuperLU_DIST 2.0), so there is no system BLAS dependency.
# ---------------------------------------------------------------------------

set(SUPERLU_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(SUPERLU_CBLAS "${CMAKE_CURRENT_SOURCE_DIR}/cblas")

# ---------------------------------------------------------------------------
# SuperLU solver sources — double-precision (d*) + shared utilities.
#
# Naming convention in SuperLU 3.0:
#   d*.c       double-precision kernels          (KEEP)
#   s*.c       single-precision kernels          (EXCLUDE)
#              — except sp_* (sparse utils) and superlu_* (timer)
#   c*.c       single-complex kernels            (EXCLUDE)
#              — except colamd.c (shared COLAMD ordering)
#   z*.c       double-complex kernels            (EXCLUDE)
#   sp_*.c     shared sparse utilities           (KEEP: sp_ienv, sp_coletree, sp_preorder)
#   superlu_*  shared timer                      (KEEP: superlu_timer.c)
#   colamd.c   shared COLAMD ordering library    (KEEP)
#   remaining  shared utilities                  (KEEP: memory, util, mmd, lsame, …)
# ---------------------------------------------------------------------------
file(GLOB ALL_SLU_SRCS "${SUPERLU_SRC}/*.c")
set(SUPERLU_SRCS "")
foreach(f ${ALL_SLU_SRCS})
    get_filename_component(fname "${f}" NAME)
    # Exclude single-precision (s*) but keep sp_* (sparse utils) and superlu_* (timer)
    if(fname MATCHES "^s" AND NOT fname MATCHES "^sp_" AND NOT fname MATCHES "^superlu_")
        continue()
    endif()
    # Exclude single-complex (c*) but keep colamd.c (shared COLAMD ordering)
    if(fname MATCHES "^c" AND NOT fname MATCHES "^colamd")
        continue()
    endif()
    # Exclude double-complex (z*)
    if(fname MATCHES "^z")
        continue()
    endif()
    # Harwell-Boeing test reader — not used by the Tahoe wrapper
    if(fname STREQUAL "dreadhb.c")
        continue()
    endif()
    list(APPEND SUPERLU_SRCS "${f}")
endforeach()

# ---------------------------------------------------------------------------
# CBLAS sources — double-precision Fortran BLAS translated by f2c.
# Provides dtrsv_, dgemm_, dgemv_, dtrsm_, daxpy_, dcopy_, idamax_, etc.
# f2c.h (types: integer, doublereal) is found via tahoe_f2c PUBLIC include dirs.
# ---------------------------------------------------------------------------
file(GLOB CBLAS_SRCS "${SUPERLU_CBLAS}/*.c")

# ---------------------------------------------------------------------------
# Single combined static library
# ---------------------------------------------------------------------------
add_library(tahoe_superlu STATIC ${SUPERLU_SRCS} ${CBLAS_SRCS})

target_include_directories(tahoe_superlu
    PUBLIC
        "${SUPERLU_SRC}"    # dsp_defs.h, supermatrix.h, util.h, Cnames.h
    PRIVATE
        "${SUPERLU_CBLAS}"  # cblas files also use Cnames.h from SUPERLU_SRC
)

# f2c.h (types: integer, doublereal, etc.) needed by CBLAS; link PRIVATE so
# the f2c runtime itself is not exposed to consumers of tahoe_superlu.
target_link_libraries(tahoe_superlu PRIVATE tahoe_f2c)

# Suppress warnings from 2005-era C and f2c-generated code
target_compile_options(tahoe_superlu PRIVATE
    -Wno-implicit-function-declaration
    -Wno-implicit-int
    -Wno-int-conversion
    -Wno-unused-result
    -Wno-return-type
    -Wno-incompatible-pointer-types
)

message(STATUS "Tahoe: SuperLU 3.0 will be built from source (bundled CBLAS, no system BLAS needed)")
