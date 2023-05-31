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

target_include_directories(${namel}_${buildl}
    PRIVATE
    ${ATHENAPP_INCLUDE_DIR}
    ${MODULES_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}
    )

set(${nameu}_LIBRARY_${buildu} ${namel}_${buildl}
    CACHE STRING
    "${namel} library ${buildl}")
mark_as_advanced(${nameu}_LIBRARY_${buildu})
