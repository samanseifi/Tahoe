# ---------------------------------------------------------------------------
# f2c runtime library (built from bundled source)
# Provides Fortran-to-C runtime routines needed by ABAQUS UMAT interfaces.
# ---------------------------------------------------------------------------

set(F2C_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(GLOB_RECURSE F2C_SRCS "${F2C_SRC_DIR}/*.c")

# Exclude:
#   UNUSED/ — platform-specific stubs not needed on Linux
#   IO/     — uses internal glibc FILE struct fields (_ptr, _flag, _base) that
#              are not accessible on glibc 2.29+. ABAQUS UMATs are pure math
#              routines and never call Fortran I/O.
list(FILTER F2C_SRCS EXCLUDE REGEX ".*/UNUSED/.*\\.c$")
list(FILTER F2C_SRCS EXCLUDE REGEX ".*/IO/.*\\.c$")

add_library(tahoe_f2c STATIC ${F2C_SRCS})

tahoe_collect_header_dirs("${F2C_SRC_DIR}" F2C_INC_DIRS)
target_include_directories(tahoe_f2c PUBLIC ${F2C_INC_DIRS})
target_compile_definitions(tahoe_f2c PUBLIC "__F2C__")

target_compile_options(tahoe_f2c PRIVATE
    -Wno-implicit-function-declaration
    -Wno-implicit-int
    -Wno-int-conversion
    -Wno-incompatible-pointer-types
    -Wno-return-type
    -Wno-unused-result
)

message(STATUS "Tahoe: f2c will be built from source (${F2C_SRC_DIR})")
