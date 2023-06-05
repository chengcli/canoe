# Installs inversion library
#
#   library: libinversion_debug.o
#   library: libinversion_release.o
#
# Define the following variables
#
#   INVERSION_INCLUDE_DIR
#   INVERSION_LIBRARY_DEBUG
#   INVERSION_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${INVERSION_INCLUDE_DIR})
#   target_link_libraries( ${INVERSION_LIBRARY_DEBUG})

set(namel inversion)
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
