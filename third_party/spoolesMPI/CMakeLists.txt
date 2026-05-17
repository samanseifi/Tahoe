# ---------------------------------------------------------------------------
# SPOOLES-MPI — MPI-parallel extension to the bundled SPOOLES 2.2 library
#
# Source: spoolesMPI/src/ from NguyenLabJHU/Tahoe-FEM (identical to the
# spoolesMPI directory shipped with SPOOLES 2.2 from Rice University).
#
# Adds distributed-memory parallel LU factorisation and solve on top of
# the serial SPOOLES data structures.  Requires MPI; works with the system
# OpenMPI selected by find_package(MPI) in the root CMakeLists.txt.
# ---------------------------------------------------------------------------

set(SPOOLESMPI_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Collect all MPI-layer C files
file(GLOB SPOOLESMPI_MPI_SRCS "${SPOOLESMPI_SRC}/MPI/*.c")

add_library(tahoe_spooles_mpi STATIC
    ${SPOOLESMPI_MPI_SRCS}
    ${SPOOLESMPI_SRC}/drivers/LU_MPI_driver.c
    ${SPOOLESMPI_SRC}/drivers/LU_MPI_driver_init.c
    ${SPOOLESMPI_SRC}/drivers/LU_MPI_driver_factorize.c
    ${SPOOLESMPI_SRC}/drivers/LU_MPI_driver_solve.c
    ${SPOOLESMPI_SRC}/drivers/LU_MPI_driver_free.c
)

# Own headers; serial SPOOLES headers (FrontMtx.h, ETree.h, …)
# come transitively through tahoe_spooles.
target_include_directories(tahoe_spooles_mpi PUBLIC
    ${SPOOLESMPI_SRC}
    ${SPOOLESMPI_SRC}/MPI
    ${SPOOLESMPI_SRC}/drivers
)

# Pull in serial SPOOLES symbols + headers, and system MPI (C interface).
# MPI::MPI_C is the C interface — SPOOLES MPI code is plain C.
target_link_libraries(tahoe_spooles_mpi PUBLIC
    tahoe_spooles
    MPI::MPI_C
)

# Expose the __SPOOLES_MPI__ define to consumers (libtahoe, SPOOLESMatrixT_mpi)
target_compile_definitions(tahoe_spooles_mpi PUBLIC __SPOOLES_MPI__)

# Silence old-C warnings present throughout SPOOLES C source base
target_compile_options(tahoe_spooles_mpi PRIVATE
    -Wno-implicit-function-declaration
    -Wno-incompatible-pointer-types
    -Wno-int-conversion
    -Wno-unused-result
    -Wno-unused-variable
    -Wno-maybe-uninitialized
)

message(STATUS "Tahoe: SPOOLES-MPI will be built from source (${SPOOLESMPI_SRC})")
