# Installs climath library
#
#   library: libclimath_debug.o
#   library: libclimath_release.o
#
# Define the following variables
#
#   CLIMATH_INCLUDE_DIR
#   CLIMATH_LIBRARY_DEBUG
#   CLIMATH_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${CLIMATH_INCLUDE_DIR})
#   target_link_libraries( ${CLIMATH_LIBRARY_DEBUG})

set(namel climath)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
  *.c
  )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
    STATIC
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_C_FLAGS_${buildu}}
    )

list(APPEND CANOE_LIBRARY_${buildu} "${namel}_${buildl}")
message(STATUS "LIBRARY = ${CANOE_LIBRARY_${buildu}}")
#set(${nameu}_LIBRARY_${buildu} ${namel}_${buildl}
#    CACHE STRING
#    "${namel} library ${buildl}")
#mark_as_advanced(${nameu}_LIBRARY_${buildu})
