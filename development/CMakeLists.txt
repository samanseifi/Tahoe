# ---------------------------------------------------------------------------
# development — research/experimental element module
#
# Compiled as a static library (libdevelopment) and linked into the tahoe
# executable via --whole-archive so that element factory registrations
# (global constructors) are preserved.
#
# The development module includes advanced element types:
#   dielectric elastomers, micromorphic, mixture, fiber composite,
#   multiscale, gradient plasticity, bio-materials, and more.
# ---------------------------------------------------------------------------

set(DEV_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Collect all .cpp files, excluding:
#   - standalone main.cpp files inside DEM_ellip3d (separate simulation codes)
#   - membrane_fluid_interaction/main.cpp
#   - *_old.cpp legacy files
#   - DEM_coupling: guarded by DEM_COUPLING_DEV which is disabled in DevelopmentElementsConfig.h;
#     its sources call FBC_CardT::AddValues/ClearValues that only exist under that guard.
file(GLOB_RECURSE DEV_SRCS "${DEV_SRC_DIR}/*.cpp")
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/DEM_ellip3d/.*\\.cpp$")
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/membrane_fluid_interaction/.*\\.cpp$")
list(FILTER DEV_SRCS EXCLUDE REGEX ".*_old\\.cpp$")
# *_continuum.cpp: legacy API variant files with mismatched signatures (treated like *_old.cpp)
list(FILTER DEV_SRCS EXCLUDE REGEX ".*_continuum\\.cpp$")
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/DEM_coupling/.*\\.cpp$")
# PMLElement: PML_ELEMENT_DEV disabled in config; abstract class issue with PMLMatT
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/PMLElement/.*\\.cpp$")
# craig_enhanced_strain_loc: ENHANCED_STRAIN_LOC_DEV_CRAIG disabled; mutual header dependency
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/craig_enhanced_strain_loc/.*\\.cpp$")
# fiber_composite: FIBER_COMP_DEV disabled in config; incomplete type errors
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/fiber_composite/.*\\.cpp$")
# meshfree_grad_plast: MESHFREE_GRAD_PLAST_DEV disabled; missing D3MeshfreeShapeFunctionT.h
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/meshfree_grad_plast/.*\\.cpp$")
# micromorphic: MICROMORPHIC_DEV disabled; combined_plasticity.cpp is a non-standalone fragment
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/micromorphic/.*\\.cpp$")
# micromorphic2: MICROMORPHIC2_DEV disabled; undeclared member variables
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/micromorphic2/.*\\.cpp$")
# micromorphic_curr_config: MICROMORPHIC_CURR_CONFIG_DEV disabled
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/micromorphic_curr_config/.*\\.cpp$")
# optimization: SOLID_OPTIMIZATION_DEV disabled; depends on fiber_composite headers
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/optimization/.*\\.cpp$")
# solid_fluid_mix: SOLID_FLUID_MIX_DEV disabled; stale variant files with undeclared members
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/solid_fluid_mix/.*\\.cpp$")
# surface_CB*: all SURFACE_CB_*_DEV disabled; missing includes / broken headers
list(FILTER DEV_SRCS EXCLUDE REGEX ".*/surface_CB[^/]*/.*\\.cpp$")

tahoe_collect_header_dirs("${DEV_SRC_DIR}" DEV_INC_DIRS)

add_library(libdevelopment STATIC ${DEV_SRCS})
set_target_properties(libdevelopment PROPERTIES OUTPUT_NAME "development")

target_include_directories(libdevelopment PUBLIC
    ${DEV_INC_DIRS}
    "${CMAKE_CURRENT_SOURCE_DIR}/config"
)

# Development depends on libtahoe (and transitively toolbox)
target_link_libraries(libdevelopment PUBLIC libtahoe)

target_compile_options(libdevelopment PRIVATE -O2)

# ---------------------------------------------------------------------------
# Inject development elements into the tahoe executable.
# We append to the existing link options rather than replacing them.
# ---------------------------------------------------------------------------
target_link_options(tahoe PRIVATE
    -Wl,--whole-archive
    $<TARGET_FILE:libdevelopment>
    -Wl,--no-whole-archive
)
target_link_libraries(tahoe PRIVATE libdevelopment)
