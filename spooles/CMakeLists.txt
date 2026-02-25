# ---------------------------------------------------------------------------
# SPOOLES 2.2 — sparse direct solver (built from bundled source)
# ---------------------------------------------------------------------------

set(SPOOLES_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(GLOB_RECURSE SPOOLES_SRCS "${SPOOLES_SRC_DIR}/*.c")

# Exclude any test programs if present (none currently in this bundled copy)
list(FILTER SPOOLES_SRCS EXCLUDE REGEX ".*/[Tt]est[s]?/.*\\.c$")

add_library(tahoe_spooles STATIC ${SPOOLES_SRCS})

tahoe_collect_header_dirs("${SPOOLES_SRC_DIR}" SPOOLES_INC_DIRS)
target_include_directories(tahoe_spooles PUBLIC ${SPOOLES_INC_DIRS})
target_compile_definitions(tahoe_spooles PUBLIC "__SPOOLES__")

# SPOOLES uses POSIX threads internally for some routines; suppress warnings
target_compile_options(tahoe_spooles PRIVATE
    -Wno-implicit-function-declaration
    -Wno-incompatible-pointer-types
    -Wno-int-conversion
    -Wno-unused-result
)

message(STATUS "Tahoe: SPOOLES will be built from source (${SPOOLES_SRC_DIR})")
