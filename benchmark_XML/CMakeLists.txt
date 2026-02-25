# ---------------------------------------------------------------------------
# benchmark_XML — regression test comparator
#   compare: reads two Tahoe output files and diffs them to tolerance
#   Depends on: toolbox only
# ---------------------------------------------------------------------------

set(COMPARATOR_SRC "${CMAKE_CURRENT_SOURCE_DIR}/comparator/src")
file(GLOB_RECURSE COMPARATOR_SRCS "${COMPARATOR_SRC}/*.cpp")
tahoe_collect_header_dirs("${COMPARATOR_SRC}" COMPARATOR_INC_DIRS)

add_executable(compare ${COMPARATOR_SRCS})
target_include_directories(compare PRIVATE ${COMPARATOR_INC_DIRS})
target_link_libraries(compare PRIVATE toolbox m)
