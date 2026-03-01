# ---------------------------------------------------------------------------
# toolbox — Tahoe utility library
#
# Provides: linear algebra containers, mesh I/O, parameter trees, geometry,
#           graph algorithms, parallel communication, and I/O streams.
#
# All Tahoe source files use flat quoted includes ("SomeName.h"), relying on
# a collected include directory.  We replicate that by adding every
# subdirectory that contains .h files to the PUBLIC include interface.
# ---------------------------------------------------------------------------

set(TOOLBOX_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Collect every .cpp file under toolbox/src/ recursively
file(GLOB_RECURSE TOOLBOX_SRCS "${TOOLBOX_SRC_DIR}/*.cpp")

# Collect all header directories (replicates the flat inc/ symlink tree)
tahoe_collect_header_dirs("${TOOLBOX_SRC_DIR}" TOOLBOX_INC_DIRS)

add_library(toolbox STATIC ${TOOLBOX_SRCS})

# PUBLIC: downstream targets (tahoe, development, ...) inherit these paths
target_include_directories(toolbox PUBLIC ${TOOLBOX_INC_DIRS})

# ---------------------------------------------------------------------------
# Compile definitions — propagated to anything that links toolbox
# ---------------------------------------------------------------------------
if(TAHOE_EXPAT)
    target_link_libraries(toolbox PUBLIC tahoe_expat)
    # __EXPAT__ already provided by the IMPORTED target's INTERFACE
endif()

if(TAHOE_SPOOLES)
    # SPOOLES headers needed at compile time for some solver wrappers
    target_compile_definitions(toolbox PUBLIC __SPOOLES__)
    target_include_directories(toolbox PUBLIC "${TAHOE_ROOT}/spooles/inc")
endif()

if(TAHOE_F2C)
    target_compile_definitions(toolbox PUBLIC __F2C__)
    target_include_directories(toolbox PUBLIC "${TAHOE_ROOT}/f2c/inc")
endif()

if(TAHOE_MPI)
    target_link_libraries(toolbox PUBLIC MPI::MPI_CXX)
    target_compile_definitions(toolbox PUBLIC __TAHOE_MPI__)
endif()

if(TAHOE_METIS)
    target_compile_definitions(toolbox PUBLIC __METIS__)
    target_include_directories(toolbox PUBLIC "${METIS_INCLUDE_DIR}")
    target_link_libraries(toolbox PUBLIC "${METIS_LIBRARY}")
endif()

if(TAHOE_SEACAS)
    target_link_libraries(toolbox PUBLIC tahoe_seacas)
endif()

# ---------------------------------------------------------------------------
# Optimisation flags
# ---------------------------------------------------------------------------
target_compile_options(toolbox PRIVATE -O2)
