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
    eos/apply_vapor_limiter.cpp
		eos/check_hydro_variables.cpp
		eos/eos.cpp
		eos/eos_hydro_ideal.cpp
		eos/${EQUATION_OF_STATE}_hydro.cpp
    #thermodynamics/*.cpp
    decomposition/*.cpp
    riemann/${RSOLVER}.cpp
    riemann/riemann.cpp
    riemann/rs_hydro_lmars.cpp
    hydro/*.cpp
    implicit/*.cpp
    reconstruct/*.cpp
    )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
    OBJECT
    ${src_files}
    )

target_include_directories(${namel}_${buildl}
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )
