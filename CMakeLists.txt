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
option(TAHOE_SEACAS   "Enable SEACAS/ExodusII mesh I/O (requires ACCESS)"  OFF)
option(TAHOE_DEV      "Build the research/development element module"      ON)

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
    add_subdirectory(expat)
endif()

if(TAHOE_SPOOLES)
    add_subdirectory(spooles)
endif()

if(TAHOE_F2C)
    add_subdirectory(f2c)
endif()

# ---------------------------------------------------------------------------
# Optional: MPI
# ---------------------------------------------------------------------------
if(TAHOE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    message(STATUS "Tahoe: MPI enabled  (${MPI_CXX_COMPILER})")
endif()

# ---------------------------------------------------------------------------
# Optional: SEACAS / ACCESS (ExodusII + NetCDF mesh I/O)
# ---------------------------------------------------------------------------
if(TAHOE_SEACAS)
    set(ACCESS_PATH "$ENV{ACCESS}" CACHE PATH "Root of SEACAS/ACCESS installation")
    if(NOT EXISTS "${ACCESS_PATH}")
        message(FATAL_ERROR
            "TAHOE_SEACAS=ON but ACCESS path '${ACCESS_PATH}' does not exist.\n"
            "Set the ACCESS environment variable or cmake -DACCESS_PATH=... or "
            "turn off TAHOE_SEACAS.")
    endif()
    add_library(tahoe_seacas INTERFACE IMPORTED)
    target_include_directories(tahoe_seacas INTERFACE "${ACCESS_PATH}/inc")
    target_link_directories   (tahoe_seacas INTERFACE "${ACCESS_PATH}/lib")
    target_link_libraries     (tahoe_seacas INTERFACE exoIIv2c netcdf)
    target_compile_definitions(tahoe_seacas INTERFACE "__ACCESS__")
    message(STATUS "Tahoe: SEACAS/ACCESS enabled  (${ACCESS_PATH})")
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
