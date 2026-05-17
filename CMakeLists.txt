cmake_minimum_required(VERSION 3.16)
project(Tahoe CXX C)

# ---------------------------------------------------------------------------
# C++ standard — use C++14 for compatibility with the old codebase
# (C++17 removes 'register' and some other C++98/03 constructs)
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Feature options
#   All optional external libraries are guarded by #ifdef in the source, so
#   disabling one simply produces a slightly reduced feature set.
# ---------------------------------------------------------------------------
option(TAHOE_EXPAT    "Use bundled expat XML parser (highly recommended)"  ON)
option(TAHOE_SPOOLES  "Use bundled SPOOLES direct sparse solver"           ON)
option(TAHOE_F2C      "Use bundled F2C (Fortran-to-C for ABAQUS UMATs)"   ON)
option(TAHOE_MPI      "Enable MPI parallelization (requires system MPI)"   OFF)
option(TAHOE_SEACAS   "Enable SEACAS/ExodusII mesh I/O (requires ACCESS)"  ON)
option(TAHOE_DEV      "Build the research/development element module"      ON)
option(TAHOE_OPENMP   "Enable OpenMP threading for explicit elements"       ON)

# ---------------------------------------------------------------------------
# OpenMP support (for ExplicitElementT batch parallelism)
# ---------------------------------------------------------------------------
if(TAHOE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "Tahoe: OpenMP enabled (${OpenMP_CXX_FLAGS})")
    else()
        message(STATUS "Tahoe: OpenMP requested but not found — disabled")
        set(TAHOE_OPENMP OFF)
    endif()
endif()

# ---------------------------------------------------------------------------
# Suppress harmless warnings from old C++98/03 code
# ---------------------------------------------------------------------------
add_compile_options(
    -Wno-deprecated
    -Wno-deprecated-declarations
    -Wno-unused-parameter
    -Wno-unused-variable
    -Wno-sign-compare
    # C++-only flags: guard with generator expression so they are not passed to
    # the C compiler when building f2c or other C-only libraries.
    $<$<COMPILE_LANGUAGE:CXX>:-Wno-reorder>
    $<$<COMPILE_LANGUAGE:CXX>:-fpermissive>
)

# ---------------------------------------------------------------------------
# Helper: collect every directory under ROOT that contains at least one
# .h file.  This replicates the flat "inc/" symlink directory that the
# old Makefile/MirrorTree.pl system used to create.
# ---------------------------------------------------------------------------
function(tahoe_collect_header_dirs ROOT RESULT_VAR)
    file(GLOB_RECURSE _headers LIST_DIRECTORIES FALSE "${ROOT}/*.h")
    set(_dirs "")
    foreach(_h ${_headers})
        get_filename_component(_d "${_h}" DIRECTORY)
        list(APPEND _dirs "${_d}")
    endforeach()
    list(REMOVE_DUPLICATES _dirs)
    set(${RESULT_VAR} "${_dirs}" PARENT_SCOPE)
endfunction()

set(TAHOE_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# ---------------------------------------------------------------------------
# Bundled libraries — built from source
#   All three are pure C and have no external dependencies.
#   Building from source removes the dependency on the old pre-built x86_64
#   .a files and fixes the PIE relocation issues they caused.
# ---------------------------------------------------------------------------

if(TAHOE_EXPAT)
    add_subdirectory(third_party/expat)
endif()

if(TAHOE_SPOOLES)
    add_subdirectory(third_party/spooles)
endif()

# ---------------------------------------------------------------------------
# Optional: SPOOLES multithreaded (single-node parallel factorisation)
#
# Adds a POSIX-threads parallel LU driver on top of the serial SPOOLES data
# structures.  No MPI required — threads only.  Requires TAHOE_SPOOLES=ON.
#
# The serial SPOOLES library also needs __POSIX_THREADS__ so that its Lock
# struct layout (pthread_mutex_t* vs void*) and Lock_init/lock/unlock match
# what the MT code expects.  This is safe for the serial solver path because
# the serial LU driver always calls SubMtxManager_init with NO_LOCK (0), so
# Lock_init never actually allocates a mutex in the serial path.
# ---------------------------------------------------------------------------
option(TAHOE_SPOOLES_MT "Build SPOOLES multithreaded solver (requires TAHOE_SPOOLES)" OFF)

if(TAHOE_SPOOLES_MT)
    if(NOT TAHOE_SPOOLES)
        message(FATAL_ERROR "TAHOE_SPOOLES_MT=ON requires TAHOE_SPOOLES=ON")
    endif()
    find_package(Threads REQUIRED)   # resolves to system -lpthread, not conda
    target_compile_definitions(tahoe_spooles PUBLIC __POSIX_THREADS__)
    add_subdirectory(third_party/spoolesMT)
    message(STATUS "Tahoe: SPOOLES-MT enabled  (${CMAKE_THREAD_LIBS_INIT})")
endif()

if(TAHOE_F2C)
    add_subdirectory(third_party/f2c)
endif()

# ---------------------------------------------------------------------------
# Optional: SuperLU 3.0 (serial sparse direct solver)
#
# Bundled SuperLU 3.0 from the original 2005 distribution. BLAS routines are
# provided by the bundled SuperLU_CBLAS (f2c-generated C), so no system BLAS
# or external dependency is required.  Requires TAHOE_F2C=ON (f2c.h types).
# ---------------------------------------------------------------------------
option(TAHOE_SUPERLU "Build bundled SuperLU 3.0 serial sparse solver" OFF)

if(TAHOE_SUPERLU)
    if(NOT TAHOE_F2C)
        message(FATAL_ERROR "TAHOE_SUPERLU=ON requires TAHOE_F2C=ON (for f2c.h types in bundled CBLAS)")
    endif()
    add_subdirectory(third_party/superlu)
    message(STATUS "Tahoe: SuperLU 3.0 enabled (bundled, no system BLAS needed)")
endif()

# ---------------------------------------------------------------------------
# Optional: MUMPS  (requires system libmumps-dev — apt-get install libmumps-dev)
#
# Uses comm_fortran=-987654 (sequential dummy communicator), so MUMPS runs
# without actual MPI communication, but the library itself still depends on
# MPI and ScaLAPACK at link time.
# ---------------------------------------------------------------------------
option(TAHOE_MUMPS "Link system MUMPS direct solver (-DTAHOE_MUMPS=ON)" OFF)
if(TAHOE_MUMPS)
    message(STATUS "Tahoe: MUMPS enabled (system library — requires libmumps-dev)")
endif()

# ---------------------------------------------------------------------------
# Optional: METIS  (requires system libmetis-dev — apt-get install libmetis-dev)
#
# Enables high-quality k-way graph partitioning for MPI domain decomposition.
# Without METIS, Tahoe uses a built-in recursive bisection (-decomp_method -0).
# With METIS, the decomposer minimises edge cuts / communication volume for
# better load balance — used automatically when __METIS__ is defined (i.e. when
# METIS is available, it is the default; disable with -no_metis at runtime).
# ---------------------------------------------------------------------------
option(TAHOE_METIS "Use system METIS graph partitioner (-DTAHOE_METIS=ON)" OFF)
if(TAHOE_METIS)
    find_path(METIS_INCLUDE_DIR metis.h
        HINTS /usr/include /usr/local/include
        DOC "METIS header directory")
    find_library(METIS_LIBRARY
        NAMES metis
        HINTS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu
        DOC "METIS library")
    if(NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
        message(FATAL_ERROR
            "TAHOE_METIS=ON but METIS not found.\n"
            "Install with: sudo apt-get install libmetis-dev\n"
            "Or set -DMETIS_INCLUDE_DIR=... -DMETIS_LIBRARY=... manually.")
    endif()
    message(STATUS "Tahoe: METIS enabled  (${METIS_LIBRARY})")
endif()

# ---------------------------------------------------------------------------
# Optional: MPI  (requires system OpenMPI/MPICH — libopenmpi-dev / libmpich-dev)
#
# Portability note: conda installs its own MPI wrappers (mpicxx, mpicc) into
# its bin/ which may appear first in PATH even when conda is not activated.
# To avoid picking up the wrong MPI, we search /usr/bin and /usr/local/bin
# explicitly before falling back to PATH.  The user can always override with
#   cmake -DMPI_CXX_COMPILER=/path/to/mpicxx
# ---------------------------------------------------------------------------
if(TAHOE_MPI)
    # Search standard system locations FIRST (NO_DEFAULT_PATH bypasses PATH,
    # so conda wrappers in PATH are ignored).  We FORCE-set the cache so that
    # stale conda values from a previous configure run are overridden — unless
    # the user has explicitly passed -DMPI_CXX_COMPILER=... on the command line
    # (CMAKE_CACHE_ARGS or -D flags always take priority over FORCE inside the
    # script because the user value is written before project() runs).
    if(NOT MPI_CXX_COMPILER OR MPI_CXX_COMPILER MATCHES "conda|miniforge|anaconda")
        find_program(_sys_mpicxx
            NAMES mpicxx mpiCC mpi-cxx
            HINTS /usr/bin /usr/local/bin /opt/local/bin
            NO_DEFAULT_PATH
            DOC "System MPI C++ wrapper (bypasses conda PATH entries)")
        if(_sys_mpicxx)
            set(MPI_CXX_COMPILER "${_sys_mpicxx}" CACHE FILEPATH
                "MPI C++ compiler wrapper (auto-detected; override with -DMPI_CXX_COMPILER=...)" FORCE)
        endif()
        unset(_sys_mpicxx CACHE)
    endif()
    if(NOT MPI_C_COMPILER OR MPI_C_COMPILER MATCHES "conda|miniforge|anaconda")
        find_program(_sys_mpicc
            NAMES mpicc mpi-cc
            HINTS /usr/bin /usr/local/bin /opt/local/bin
            NO_DEFAULT_PATH)
        if(_sys_mpicc)
            set(MPI_C_COMPILER "${_sys_mpicc}" CACHE FILEPATH
                "MPI C compiler wrapper (auto-detected)" FORCE)
        endif()
        unset(_sys_mpicc CACHE)
    endif()

    find_package(MPI REQUIRED COMPONENTS CXX C)

    # Warn if a non-system wrapper slipped through (e.g. conda not deactivated)
    if(MPI_CXX_COMPILER MATCHES "conda|miniforge|anaconda")
        message(WARNING
            "MPI wrapper looks like a conda install: ${MPI_CXX_COMPILER}\n"
            "This may conflict with system SEACAS/HDF5 libraries.\n"
            "Override: cmake -DMPI_CXX_COMPILER=/usr/bin/mpicxx")
    endif()

    message(STATUS "Tahoe: MPI enabled  (${MPI_CXX_COMPILER})")
    add_subdirectory(third_party/spoolesMPI)
endif()

# ---------------------------------------------------------------------------
# Optional: SEACAS / ACCESS (ExodusII + NetCDF mesh I/O)
#
# Two discovery modes (tried in order):
#   1. SEACAS/ACCESS tree  — set ACCESS_PATH or the ACCESS env var
#   2. System packages     — headers/libs found via find_path / find_library
#      (e.g. apt install libexodusii-dev libnetcdf-dev on Ubuntu/Debian)
# ---------------------------------------------------------------------------
if(TAHOE_SEACAS)
    add_library(tahoe_seacas INTERFACE IMPORTED)
    target_compile_definitions(tahoe_seacas INTERFACE "__ACCESS__")

    set(ACCESS_PATH "$ENV{ACCESS}" CACHE PATH
        "Root of SEACAS/ACCESS installation (optional; falls back to system search)")

    if(EXISTS "${ACCESS_PATH}")
        # --- Mode 1: SEACAS/ACCESS tree ---
        target_include_directories(tahoe_seacas INTERFACE "${ACCESS_PATH}/inc")
        target_link_directories   (tahoe_seacas INTERFACE "${ACCESS_PATH}/lib")
        target_link_libraries     (tahoe_seacas INTERFACE exoIIv2c netcdf)
        message(STATUS "Tahoe: SEACAS/ACCESS enabled via ACCESS tree (${ACCESS_PATH})")
    else()
        # --- Mode 2: System-installed packages ---
        find_path(EXODUSII_INCLUDE_DIR  NAMES exodusII.h
                  HINTS /usr/include /usr/local/include)
        find_library(EXODUSII_LIB       NAMES exoIIv2c ExodusII exodus
                     HINTS /usr/lib /usr/local/lib
                           /usr/lib/x86_64-linux-gnu
                           /usr/lib/aarch64-linux-gnu)
        find_library(NETCDF_LIB         NAMES netcdf
                     HINTS /usr/lib /usr/local/lib
                           /usr/lib/x86_64-linux-gnu
                           /usr/lib/aarch64-linux-gnu)

        if(NOT EXODUSII_INCLUDE_DIR OR NOT EXODUSII_LIB OR NOT NETCDF_LIB)
            message(FATAL_ERROR
                "TAHOE_SEACAS=ON but neither an ACCESS tree nor system ExodusII "
                "packages were found.\n"
                "  • Set ACCESS_PATH / ACCESS env var, or\n"
                "  • Install system packages (e.g. sudo apt install libexodusii-dev libnetcdf-dev), or\n"
                "  • Disable: cmake -DTAHOE_SEACAS=OFF")
        endif()

        target_include_directories(tahoe_seacas INTERFACE "${EXODUSII_INCLUDE_DIR}")
        target_link_libraries     (tahoe_seacas INTERFACE "${EXODUSII_LIB}" "${NETCDF_LIB}")
        message(STATUS "Tahoe: SEACAS/ACCESS enabled via system packages")
        message(STATUS "  ExodusII header : ${EXODUSII_INCLUDE_DIR}")
        message(STATUS "  ExodusII lib    : ${EXODUSII_LIB}")
        message(STATUS "  NetCDF lib      : ${NETCDF_LIB}")
    endif()
endif()

# ---------------------------------------------------------------------------
# Sub-modules
# ---------------------------------------------------------------------------
add_subdirectory(toolbox)
add_subdirectory(tahoe)         # the 'tahoe' sub-directory (the FEM module)

if(TAHOE_DEV)
    add_subdirectory(development)
endif()

# ---------------------------------------------------------------------------
# Contrib utilities and benchmark comparator
# ---------------------------------------------------------------------------
add_subdirectory(contrib)
add_subdirectory(benchmark_XML)

# ---------------------------------------------------------------------------
# Unit tests (Google Test, fetched automatically)
# ---------------------------------------------------------------------------
option(TAHOE_TESTS "Build Google Test unit tests" ON)
if(TAHOE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
