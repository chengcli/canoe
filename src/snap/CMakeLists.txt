# Installs snap library
#
#   library: libsnap_debug.o
#   library: libsnap_release.o
#
# Define the following variables
#
#   SNAP_INCLUDE_DIR
#   SNAP_LIBRARY_DEBUG
#   SNAP_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${SNAP_INCLUDE_DIR})
#   target_link_libraries( ${SNAP_LIBRARY_DEBUG})

set(namel snap)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
    *.cpp
    eos/*.cpp
    thermodynamics/*.cpp
    decomposition/*.cpp
    reconstruct/*.cpp
    implicit/*.cpp
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

set(${nameu}_LIBRARY_${buildu} ${namel}_${buildl}
    CACHE STRING
    "${namel} library ${buildl}")
mark_as_advanced(${nameu}_LIBRARY_${buildu})
