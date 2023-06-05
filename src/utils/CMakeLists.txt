# Installs utils library
#
#   library: libutils_debug.o
#   library: libutils_release.o
#
# Define the following variables
#
#   UTILS_INCLUDE_DIR
#   UTILS_LIBRARY_DEBUG
#   UTILS_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${UTILS_INCLUDE_DIR})
#   target_link_libraries( ${UTILS_LIBRARY_DEBUG})

set(namel utils)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
    *.cpp
    )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
    OBJECT
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )
