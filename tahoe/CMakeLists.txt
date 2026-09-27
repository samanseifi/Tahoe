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

# ---------------------------------------------------------------------------
# MUMPS (default linear solver) and the optimized BLAS it runs on (#80).
# Neither is required: a missing library gives a warning and a fallback.
# ---------------------------------------------------------------------------
set(TAHOE_HAVE_MUMPS OFF)
set(TAHOE_BLAS_LIBRARIES "")
if(TAHOE_MUMPS)
    # System MUMPS carries its own transitive deps (ScaLAPACK, MPI Fortran,
    # gfortran); we only need the C MPI headers for MUMPSMatrixT_mpi and the
    # MPI_Init helper.  The Fortran MPI component would require enabling
    # Fortran in CMake, so it is not requested.
    find_package(MPI QUIET COMPONENTS C)

    set(_tahoe_mumps_hints
        /usr/lib /usr/lib64 /usr/local/lib /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/openmpi/lib)
    find_library(MUMPS_D_LIBRARY   NAMES dmumps       HINTS ${_tahoe_mumps_hints})
    find_library(MUMPS_CMN_LIBRARY NAMES mumps_common HINTS ${_tahoe_mumps_hints})
    find_library(PORD_LIBRARY      NAMES pord         HINTS ${_tahoe_mumps_hints})
    find_library(SCALAPACK_LIBRARY NAMES scalapack-openmpi scalapack HINTS ${_tahoe_mumps_hints})
    find_path(MUMPS_INCLUDE_DIR    NAMES dmumps_c.h
        HINTS /usr/include /usr/local/include /opt/homebrew/include
        PATH_SUFFIXES mumps MUMPS)

    if(MUMPS_D_LIBRARY AND MUMPS_CMN_LIBRARY AND MUMPS_INCLUDE_DIR AND MPI_C_FOUND)
        set(TAHOE_HAVE_MUMPS ON)
        target_include_directories(libtahoe PUBLIC ${MUMPS_INCLUDE_DIR})
        target_compile_definitions(libtahoe PUBLIC __MUMPS__)
        target_link_libraries(libtahoe PUBLIC ${MUMPS_D_LIBRARY} ${MUMPS_CMN_LIBRARY} MPI::MPI_C)
        foreach(_lib PORD_LIBRARY SCALAPACK_LIBRARY)
            if(${_lib})
                target_link_libraries(libtahoe PUBLIC ${${_lib}})
            endif()
        endforeach()
    else()
        message(WARNING
            "Tahoe: MUMPS not found (dmumps=${MUMPS_D_LIBRARY}, header=${MUMPS_INCLUDE_DIR}, "
            "MPI C=${MPI_C_FOUND}).  Building without MUMPS; SPOOLES is the default solver.\n"
            "For the faster default solver install MUMPS and OpenBLAS, e.g. on Debian/Ubuntu:\n"
            "    sudo apt-get install libmumps-dev libopenblas-dev\n"
            "or configure with -DTAHOE_MUMPS=OFF to silence this warning.")
    endif()

    if(TAHOE_HAVE_MUMPS AND TAHOE_BLAS)
        # Search system and CMAKE_PREFIX_PATH locations only: PATH-derived
        # prefixes (for example an active conda environment) would put that
        # environment's lib directory on the RPATH, where its MPI, curl and
        # gfortran shadow the system libraries MUMPS was built against.
        # Point at another BLAS with -DCMAKE_PREFIX_PATH=... or -DBLAS_LIBRARIES=...
        set(_tahoe_saved_env_path ${CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH})
        set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH OFF)
        set(BLA_VENDOR ${TAHOE_BLAS_VENDOR})
        find_package(BLAS QUIET)
        set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH ${_tahoe_saved_env_path})
        unset(BLA_VENDOR)
        if(BLAS_FOUND AND BLAS_LIBRARIES)
            set(TAHOE_BLAS_LIBRARIES ${BLAS_LIBRARIES})
        else()
            message(WARNING
                "Tahoe: ${TAHOE_BLAS_VENDOR} BLAS not found; MUMPS will use the BLAS it was built "
                "against, which on many distributions is the unoptimized reference BLAS (about 6x "
                "slower).  Install it (Debian/Ubuntu: sudo apt-get install libopenblas-dev), choose "
                "another vendor with -DTAHOE_BLAS_VENDOR=..., or set -DTAHOE_BLAS=OFF.")
        endif()
    endif()
endif()
set(TAHOE_HAVE_MUMPS ${TAHOE_HAVE_MUMPS} PARENT_SCOPE)
set(TAHOE_BLAS_LIBRARIES "${TAHOE_BLAS_LIBRARIES}" PARENT_SCOPE)

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

# Link the optimized BLAS directly into the executable (#80).  It must be a
# direct dependency (DT_NEEDED of tahoe itself), because the loader resolves
# symbols in breadth-first order: a BLAS needed by the executable wins over the
# libblas.so.3 that libdmumps pulls in one level down.  tahoe calls no BLAS
# itself, so --no-as-needed keeps the linker from dropping the library.
if(TAHOE_BLAS_LIBRARIES)
    if(APPLE)
        target_link_libraries(tahoe PRIVATE ${TAHOE_BLAS_LIBRARIES})
    else()
        target_link_libraries(tahoe PRIVATE
            "-Wl,--push-state,--no-as-needed" ${TAHOE_BLAS_LIBRARIES} "-Wl,--pop-state")
    endif()
endif()

install(TARGETS tahoe DESTINATION bin)
