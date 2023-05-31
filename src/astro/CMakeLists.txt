# Installs cliastro library
#
#   library: libcliastro_debug.o
#   library: libcliastro_release.o
#
# Define the following variables
#
#   ASTRO_INCLUDE_DIR
#   ASTRO_LIBRARY_DEBUG
#   ASTRO_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${ASTRO_INCLUDE_DIR})
#   target_link_libraries( ${ASTRO_LIBRARY_DEBUG})

set(namel astro)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
  *.cpp
  Jupiter/*.cpp
  Saturn/*.cpp
  )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
    STATIC
    ${src_files})

  #add_dependencies(cliastro_${buildl} athenapp_${buildl})

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

target_link_libraries(${namel}_${buildl}
    ${CLIMATH_LIBRARY_${buildu}}
    ${DEBUGGER_LIBRARY_${buildu}}
    ${ATHENAPP_LIBRARY_${buildu}}
    )

set(${nameu}_LIBRARY_${buildu} ${namel}_${buildl}
    CACHE STRING
    "${namel} library ${buildl}")
mark_as_advanced(${nameu}_LIBRARY_${buildu})
