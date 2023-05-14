# Installs climath library
#
#   library: libclimath_debug.o
#   library: libclimath_release.o
#
# Define the following variables
#
#   CLIMATH_FOUND
#   CLIMATH_INCLUDE_DIR
#   CLIMATH_LIBRARY_DEBUG
#   CLIMATH_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${CLIMATH_INCLUDE_DIR})
#   target_link_libraries( ${CLIMATH_LIBRARY_DEBUG})
#
cmake_minimum_required(VERSION 3.10)

project(
  climath
  LANGUAGES C
  )
set(CMAKE_C_STANDARD 99)
message(STATUS "")
message(STATUS "== Setting up climath library ==")

file(GLOB _src_files
  *.c
  )

foreach(build ${BUILD_TYPES})
  string(TOLOWER ${build} buildl)
  add_library(climath_${buildl}
    ${_src_files})

  set_target_properties(climath_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_C_FLAGS_${build}}
    )

  set(CLIMATH_LIBRARY_${build} climath_${buildl}
    CACHE STRING 
    "climath library ${build}")
  mark_as_advanced(CLIMATH_LIBRARY_${build})
endforeach()

set(CLIMATH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
  CACHE PATH 
  "climath include directory")
set(CLIMATH_FOUND TRUE CACHE BOOL "FIND CLIMATH LIBRARY")
mark_as_advanced(CLIMATH_INCLUDE_DIR)
mark_as_advanced(CLIMATH_FOUND)
