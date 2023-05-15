# Installs outputs library
#
#   library: liboutputs_debug.o
#   library: liboutputs_release.o
#
# Define the following variables
#
#   OUTPUTS_INCLUDE_DIR
#   OUTPUTS_LIBRARY_DEBUG
#   OUTPUTS_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${OUTPUTS_INCLUDE_DIR})
#   target_link_libraries( ${OUTPUTS_LIBRARY_DEBUG})

set(namel outputs)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
    *.cpp
    )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
	STATIC
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

target_include_directories(${namel}_${buildl}
	PRIVATE 
    ${ATHENAPP_INCLUDE_DIR}
	${CMAKE_BINARY_DIR}
	)

set(${nameu}_LIBRARY_${buildu} ${nameu}_${buildl}
    CACHE STRING 
    "${namel} library ${buildl}")
mark_as_advanced(${nameu}_LIBRARY_${buildu})
