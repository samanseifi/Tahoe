# ---------------------------------------------------------------------------
# tests/CMakeLists.txt — Google Test unit tests for Tahoe
# ---------------------------------------------------------------------------

include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
    GIT_SHALLOW    TRUE
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)
enable_testing()

# ---------------------------------------------------------------------------
# Helper: toolbox-only test (fast, no FEM infrastructure)
# ---------------------------------------------------------------------------
macro(tahoe_gtest target_name source_file)
    add_executable(${target_name} ${source_file})
    target_link_libraries(${target_name} PRIVATE GTest::gtest_main toolbox)
    target_compile_options(${target_name} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-fpermissive>
        -Wno-deprecated
    )
    gtest_discover_tests(${target_name})
endmacro()

# ---------------------------------------------------------------------------
# Helper: full libtahoe test (material models, needs whole-archive)
# ---------------------------------------------------------------------------
macro(tahoe_gtest_full target_name source_file)
    add_executable(${target_name} ${source_file})
    target_link_libraries(${target_name} PRIVATE GTest::gtest_main libtahoe m dl)
    target_link_options(${target_name} PRIVATE
        -Wl,--whole-archive
        $<TARGET_FILE:libtahoe>
        $<TARGET_FILE:toolbox>
        -Wl,--no-whole-archive
    )
    target_compile_options(${target_name} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-fpermissive>
        -Wno-deprecated
    )
    gtest_discover_tests(${target_name})
endmacro()

# ---------------------------------------------------------------------------
# Toolbox-only tests
# ---------------------------------------------------------------------------
tahoe_gtest(test_dArrayT     toolbox/test_dArrayT.cpp)
tahoe_gtest(test_dMatrixT    toolbox/test_dMatrixT.cpp)
tahoe_gtest(test_dSymMatrixT toolbox/test_dSymMatrixT.cpp)
tahoe_gtest(test_QuadT       geometry/test_QuadT.cpp)

# ---------------------------------------------------------------------------
# Full libtahoe tests (material models)
# ---------------------------------------------------------------------------
tahoe_gtest_full(test_IsotropicT  materials/test_IsotropicT.cpp)
tahoe_gtest_full(test_NeoHookean  materials/test_NeoHookean.cpp)
