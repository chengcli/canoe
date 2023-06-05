# Installs debugger library
#
#   library: libdebugger_debug.o
#   library: libdebugger_release.o
#
# Define the following variables
#
#   DEBUGGER_INCLUDE_DIR
#   DEBUGGER_LIBRARY_DEBUG
#   DEBUGGER_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${DEBUGGER_INCLUDE_DIR})
#   target_link_libraries( ${DEBUGGER_LIBRARY_DEBUG})

set(namel debugger)
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
