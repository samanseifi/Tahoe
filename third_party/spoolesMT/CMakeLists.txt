# ---------------------------------------------------------------------------
# SPOOLES-MT — multithreaded extension to the bundled SPOOLES 2.2 library
#
# Source: spoolesMT/src/ from NguyenLabJHU/Tahoe-FEM (identical to the
# spoolesMT directory shipped with the original SPOOLES 2.2 distribution).
#
# The library adds POSIX-thread-parallel factorisation and solve on top of
# the serial SPOOLES data structures.  It does NOT replace the serial solver;
# both can coexist in the same executable.
#
# Required defines (set here and propagated to consumers):
#   __POSIX_THREADS__  — selects the pthreads code path inside SPOOLES Lock.h
#   (tahoe_spooles itself also needs __POSIX_THREADS__ so that its Lock object
#    layout matches — see the note in the root CMakeLists.txt)
# ---------------------------------------------------------------------------

set(SPOOLESMT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(tahoe_spooles_mt STATIC
    ${SPOOLESMT_SRC}/MT/MT.factorMT.c
    ${SPOOLESMT_SRC}/MT/MT.solveMT.c
    ${SPOOLESMT_SRC}/MT/MT.mvm.c
    ${SPOOLESMT_SRC}/MT/MT.QRfactorMT.c
    ${SPOOLESMT_SRC}/MT/MT.QRsolveMT.c
    ${SPOOLESMT_SRC}/drivers/LU_MT_driver.c
    ${SPOOLESMT_SRC}/drivers/LU_MT_driver_init.c
    ${SPOOLESMT_SRC}/drivers/LU_MT_driver_factorize.c
    ${SPOOLESMT_SRC}/drivers/LU_MT_driver_solve.c
    ${SPOOLESMT_SRC}/drivers/LU_MT_driver_free.c
)

# Own headers; serial SPOOLES headers (FrontMtx.h, ETree.h, Lock.h, …)
# come transitively through tahoe_spooles.
target_include_directories(tahoe_spooles_mt PUBLIC
    ${SPOOLESMT_SRC}
    ${SPOOLESMT_SRC}/MT
    ${SPOOLESMT_SRC}/drivers
)

# Pull in serial SPOOLES symbols + headers, and system pthreads.
# Threads::Threads resolves to -pthread on Linux (system libpthread, not conda).
target_link_libraries(tahoe_spooles_mt PUBLIC
    tahoe_spooles
    Threads::Threads
)

# Select the POSIX-threads code path in Lock.h / factorMT.c / solveMT.c / mvm.c
target_compile_definitions(tahoe_spooles_mt PUBLIC __POSIX_THREADS__)

# Silence old-C warnings that exist throughout the SPOOLES C source base
target_compile_options(tahoe_spooles_mt PRIVATE
    -Wno-implicit-function-declaration
    -Wno-incompatible-pointer-types
    -Wno-int-conversion
    -Wno-unused-result
    -Wno-unused-variable
    -Wno-maybe-uninitialized
)

message(STATUS "Tahoe: SPOOLES-MT will be built from source (${SPOOLESMT_SRC})")
