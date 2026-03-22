# ---------------------------------------------------------------------------
# tahoe — main Tahoe FEM module
#
# Builds:
#   libtahoe  — static library with all elements, solvers, integrators, nodes
#   tahoe     — executable (links libtahoe + toolbox + bundled libs)
# ---------------------------------------------------------------------------

set(TAHOE_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Collect all .cpp and .c files under tahoe/src/, excluding:
#   - the entry-point main.cpp
#   - *_old.cpp files (legacy/superseded implementations kept for reference)
file(GLOB_RECURSE ALL_TAHOE_SRCS "${TAHOE_SRC_DIR}/*.cpp" "${TAHOE_SRC_DIR}/*.c")
list(FILTER ALL_TAHOE_SRCS EXCLUDE REGEX ".*/main/main\\.cpp$")
list(FILTER ALL_TAHOE_SRCS EXCLUDE REGEX ".*_old\\.cpp$")

# Collect all header directories
tahoe_collect_header_dirs("${TAHOE_SRC_DIR}" TAHOE_INC_DIRS)

# ---------------------------------------------------------------------------
# Static library:  libtahoe
# ---------------------------------------------------------------------------
add_library(libtahoe STATIC ${ALL_TAHOE_SRCS})
set_target_properties(libtahoe PROPERTIES OUTPUT_NAME "tahoe")

# PUBLIC includes: both tahoe's own headers and toolbox headers
# (toolbox headers come in automatically via target_link_libraries below)
# Also include the pre-generated config headers (ElementsConfig.h, SolidMaterialsConfig.h, etc.)
# and all development headers (some tahoe sources directly include development headers),
# replacing the broken flat inc/ symlink mirror from the old tahoe-manager build.
tahoe_collect_header_dirs("${TAHOE_ROOT}/development/src" DEV_SRC_INC_DIRS)

target_include_directories(libtahoe PUBLIC
    ${TAHOE_INC_DIRS}
    "${TAHOE_ROOT}/tahoe/config"
    "${TAHOE_ROOT}/development/config"
    ${DEV_SRC_INC_DIRS}
)

# __DEVELOPMENT__ activates the #include "DevelopmentElementsConfig.h" guards
# throughout libtahoe sources (ElementListT, FieldT, SolverT, NodeManagerT, etc.)
# so that development element flags (DIELECTRIC_ELASTOMER, etc.) resolve correctly.
target_compile_definitions(libtahoe PUBLIC __DEVELOPMENT__)

# toolbox is a PUBLIC dependency — its include dirs are transitively exported
target_link_libraries(libtahoe PUBLIC toolbox)

if(TAHOE_SPOOLES)
    target_link_libraries(libtahoe PUBLIC tahoe_spooles)
endif()

if(TAHOE_SPOOLES_MT)
    # Pull in the MT library (which already links tahoe_spooles + Threads::Threads)
    target_link_libraries(libtahoe PUBLIC tahoe_spooles_mt)
    # Activate SPOOLESMatrixT_MT in SolverT.cpp and SPOOLESMatrixT_MT.cpp
    target_compile_definitions(libtahoe PUBLIC __SPOOLES_MT__ __POSIX_THREADS__)
endif()

if(TAHOE_F2C)
    target_link_libraries(libtahoe PUBLIC tahoe_f2c)
endif()

if(TAHOE_SUPERLU)
    target_link_libraries(libtahoe PUBLIC tahoe_superlu)
    # Activates SuperLUMatrixT in SolverT.cpp and SuperLUMatrixT.cpp
    target_compile_definitions(libtahoe PUBLIC __SUPERLU__)
endif()

if(TAHOE_MPI)
    target_link_libraries(libtahoe PUBLIC MPI::MPI_CXX tahoe_spooles_mpi)
    # __SPOOLES_MPI__ comes transitively from tahoe_spooles_mpi PUBLIC defines
endif()

if(TAHOE_MUMPS)
    # MUMPS from the system package manager (e.g. libmumps-dev on Ubuntu).
    # libdmumps.so is a shared library that already carries transitive deps on
    # libscalapack, libmpi_mpifh, and libgfortran — so we only need to link
    # MPI::MPI_C for the C-level MPI headers used by MUMPSMatrixT_mpi.
    # We do NOT request the Fortran MPI component because that requires the
    # Fortran language to be enabled in CMake.
    find_package(MPI REQUIRED COMPONENTS C)

    find_library(MUMPS_D_LIBRARY   NAMES dmumps
        HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/openmpi/lib)
    find_library(MUMPS_CMN_LIBRARY NAMES mumps_common
        HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/openmpi/lib)
    find_library(PORD_LIBRARY      NAMES pord
        HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/openmpi/lib)
    find_library(SCALAPACK_LIBRARY NAMES scalapack-openmpi scalapack
        HINTS /usr/lib /usr/lib/x86_64-linux-gnu)
    find_path(MUMPS_INCLUDE_DIR    NAMES dmumps_c.h
        HINTS /usr/include /usr/include/mumps_seq)

    if(NOT MUMPS_D_LIBRARY OR NOT MUMPS_INCLUDE_DIR)
        message(FATAL_ERROR
            "TAHOE_MUMPS=ON but MUMPS libraries not found.\n"
            "Install with: sudo apt-get install libmumps-dev")
    endif()

    target_include_directories(libtahoe PUBLIC ${MUMPS_INCLUDE_DIR})
    target_compile_definitions(libtahoe PUBLIC __MUMPS__)
    target_link_libraries(libtahoe PUBLIC
        ${MUMPS_D_LIBRARY}
        ${MUMPS_CMN_LIBRARY}
        ${PORD_LIBRARY}
        ${SCALAPACK_LIBRARY}
        MPI::MPI_C
    )
endif()

if(TAHOE_SEACAS)
    target_link_libraries(libtahoe PUBLIC tahoe_seacas)
endif()

if(TAHOE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(libtahoe PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_options(libtahoe PRIVATE -O2)

# ---------------------------------------------------------------------------
# Executable: tahoe
# ---------------------------------------------------------------------------
add_executable(tahoe "${TAHOE_SRC_DIR}/main/main.cpp")

target_link_libraries(tahoe PRIVATE
    libtahoe
    # system libraries often needed by old UNIX code
    m       # math
    dl      # dynamic loading
)

# Link against the whole archive of libtahoe so that static initializers
# (element/material factories registered via global constructors) are kept.
# Without this, the linker may discard unreferenced translation units and
# element types will silently disappear at runtime.
#
target_link_options(tahoe PRIVATE
    -Wl,--whole-archive
    $<TARGET_FILE:libtahoe>
    $<TARGET_FILE:toolbox>
    -Wl,--no-whole-archive
    # Prevent SuperLU's bundled BLAS routines (idamax_, dgemm_, etc.) from
    # polluting the binary's dynamic symbol table.  Without this flag, shared
    # libraries like libdmumps resolve idamax_/dgemm_/… from our binary
    # instead of from their linked system BLAS (libblas.so.3), causing
    # crashes inside MUMPS's frontal factorization kernel.
    -Wl,--exclude-libs,ALL
)

install(TARGETS tahoe DESTINATION bin)
