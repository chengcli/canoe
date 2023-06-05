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
    OBJECT
    ${src_files})

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )
