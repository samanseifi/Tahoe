# ---------------------------------------------------------------------------
# contrib — utility executables that ship alongside Tahoe
#
# All targets land in the same CMAKE_RUNTIME_OUTPUT_DIRECTORY as tahoe
# (i.e. build/bin/).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Helper: collect headers for each contrib tool's own source tree
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# MakeCSE — mesh pre-processor that inserts cohesive surface elements
#   Depends on: toolbox + libtahoe (uses ModelManagerT, geometry classes, etc.)
# ---------------------------------------------------------------------------
set(MAKECSE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/MakeCSE/src")
file(GLOB_RECURSE MAKECSE_SRCS "${MAKECSE_SRC}/*.cpp")
tahoe_collect_header_dirs("${MAKECSE_SRC}" MAKECSE_INC_DIRS)

add_executable(MakeCSE ${MAKECSE_SRCS})
target_include_directories(MakeCSE PRIVATE ${MAKECSE_INC_DIRS})
target_link_libraries(MakeCSE PRIVATE libtahoe)
target_link_options(MakeCSE PRIVATE
    -Wl,--whole-archive
    $<TARGET_FILE:libtahoe>
    $<TARGET_FILE:toolbox>
    -Wl,--no-whole-archive
)

# ---------------------------------------------------------------------------
# translate — converts Tahoe output between formats; extracts nodal/element data
#   Depends on: toolbox + libtahoe (uses ModelManagerT, I/O managers)
# ---------------------------------------------------------------------------
set(TRANSLATE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/translate/src")
file(GLOB_RECURSE TRANSLATE_SRCS "${TRANSLATE_SRC}/*.cpp")
tahoe_collect_header_dirs("${TRANSLATE_SRC}" TRANSLATE_INC_DIRS)

add_executable(translate ${TRANSLATE_SRCS})
target_include_directories(translate PRIVATE ${TRANSLATE_INC_DIRS})
target_link_libraries(translate PRIVATE libtahoe)
target_link_options(translate PRIVATE
    -Wl,--whole-archive
    $<TARGET_FILE:libtahoe>
    $<TARGET_FILE:toolbox>
    -Wl,--no-whole-archive
)

# ---------------------------------------------------------------------------
# wrap — wraps/converts mesh files using Tahoe's I/O layer
#   Depends on: toolbox + libtahoe
# ---------------------------------------------------------------------------
set(WRAP_SRC "${CMAKE_CURRENT_SOURCE_DIR}/wrap/src")
file(GLOB_RECURSE WRAP_SRCS "${WRAP_SRC}/*.cpp")
tahoe_collect_header_dirs("${WRAP_SRC}" WRAP_INC_DIRS)

add_executable(wrap ${WRAP_SRCS})
target_include_directories(wrap PRIVATE ${WRAP_INC_DIRS})
target_link_libraries(wrap PRIVATE libtahoe)
target_link_options(wrap PRIVATE
    -Wl,--whole-archive
    $<TARGET_FILE:libtahoe>
    $<TARGET_FILE:toolbox>
    -Wl,--no-whole-archive
)

# ---------------------------------------------------------------------------
# Dakota post-processing tools — standalone C++ utilities using toolbox only
# ---------------------------------------------------------------------------
set(DAKOTA_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Dakota")

foreach(_tool area_diff calc_disp_cost extract_1D extract_CT extract_Disp
              least_square least_square_trap)
    if(EXISTS "${DAKOTA_SRC}/${_tool}.cpp")
        add_executable(${_tool} "${DAKOTA_SRC}/${_tool}.cpp")
        target_link_libraries(${_tool} PRIVATE toolbox m)
    endif()
endforeach()

# dpp is a Perl script, not compiled — skip
