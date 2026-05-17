# ---------------------------------------------------------------------------
# expat 1.95.7 — lightweight XML parser (built from bundled source)
# ---------------------------------------------------------------------------

set(EXPAT_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/expat/expat-1.95.7/lib")

add_library(tahoe_expat STATIC
    ${EXPAT_SRC_DIR}/xmlparse.c
    ${EXPAT_SRC_DIR}/xmlrole.c
    ${EXPAT_SRC_DIR}/xmltok.c
)

target_include_directories(tahoe_expat PUBLIC ${EXPAT_SRC_DIR})
target_compile_definitions(tahoe_expat PUBLIC "__EXPAT__")

# expat needs XML_NS and XML_DTD to enable full feature set (matching original build)
# HAVE_MEMMOVE tells xmlparse.c that memmove() is available (it always is on Linux)
target_compile_definitions(tahoe_expat PRIVATE XML_NS XML_DTD HAVE_MEMMOVE)

# Suppress warnings from ancient C code
target_compile_options(tahoe_expat PRIVATE
    -Wno-implicit-function-declaration
    -Wno-unused-result
)

message(STATUS "Tahoe: expat will be built from source (${EXPAT_SRC_DIR})")
