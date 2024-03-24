cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(canoe
  LANGUAGES CXX C Fortran
	)

message(STATUS "")
message(STATUS "== Setting up canoe library ==")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Default Build Type = Release")
  message(STATUS "")
else()
  message(STATUS "Build Type = ${CMAKE_BUILD_TYPE}")
  message(STATUS "")
endif()

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

## 1. set up compiler flags ##
message(STATUS "1. Set up project compiler flags ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/compilers.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/compilers.cmake)

## 2. set up project specific configuration ##
message(STATUS "2. Set up project parameters ...")

# load custom task
message(STATUS "Load custom task = ${TASK}")
if (NOT ${TASK} STREQUAL "")
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/examples/${TASK}.cmake)
    message(STATUS "Load custom setting - ${TASK} -")
    message(STATUS "Include ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${TASK}.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/examples/${TASK}.cmake)
  else()
    message(FATAL_ERROR "Cannot find custom file cmake/${TASK}.cmake")
  endif()
  set(TASK ${TASK} CACHE STRING "Custom task to be performed")
else()
  set(TASK "default" CACHE STRING "Custom task to be performed")
endif()

message(STATUS "Include ${CMAKE_CURRENT_SOURCE_DIR}/cmake/parameters.cmake")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/parameters.cmake)
find_package(Eigen3 REQUIRED)

## 3. set up project system libraries ##
message(STATUS "3. Set up system libraries")
#message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/athena.cmake")
#include(${CMAKE_SOURCE_DIR}/cmake/eigen.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/yamlpp.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/gtest.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/athena.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/application.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/rfm.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/disort.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/rrtmg.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/pvfmm.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/pybind11.cmake)

## 4. set up project configure file and library ##
message(STATUS "4. Set up project libraries")
configure_file(${CMAKE_SOURCE_DIR}/configure.hpp.in configure.hpp @ONLY)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(data)

## 5. set up python binding ##
if (PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

## 6. set up unit tests and examples
message(STATUS "6. Set up unit tests")

# load all macros
FILE(GLOB _macro_files "${CMAKE_SOURCE_DIR}/cmake/macros/*.cmake")
FOREACH(_file ${_macro_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

add_subdirectory(examples)
add_subdirectory(tests)

## 7. add uninstall targets ##
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(canoe_uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
