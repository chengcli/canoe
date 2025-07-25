cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(canoe
  LANGUAGES CXX C Fortran
	)

# search for <Package>_ROOT
cmake_policy(SET CMP0074 NEW)

find_program(Python3_EXECUTABLE NAMES python)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import torch; print(int(torch._C._GLIBCXX_USE_CXX11_ABI))"
  RESULT_VARIABLE _ABI_RESULT
  OUTPUT_VARIABLE _ABI_OUTPUT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT _ABI_RESULT EQUAL 0)
  message(FATAL_ERROR "Failed to detect Torch GLIBCXX ABI")
endif()

# Set the _GLIBCXX_USE_CXX11_ABI flag
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=${_ABI_OUTPUT})

# Set MacOS deployment target
set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0")

string(ASCII 27 Esc)
set(ColorReset "${Esc}[m")
set(Green      "${Esc}[32m")
set(Blue       "${Esc}[34m")

message(STATUS "")
message(STATUS "${Green}==== ${PROJECT_NAME} configure begin ====${ColorReset}")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Build Type = Release")
else()
  message(STATUS "Build Type = ${CMAKE_BUILD_TYPE}")
endif()

# load all modules
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

# load all macros
file(GLOB _macro_files "${CMAKE_SOURCE_DIR}/cmake/macros/*.cmake")
foreach(_file ${_macro_files})
  MESSAGE(STATUS "Include ${_file}")
  include(${_file})
endforeach()

## 1. set up compiler flags ##
message(STATUS "${Blue}${PROJECT_NAME}-1. Setting up project compiler flags ...${ColorReset}")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/compilers.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/compilers.cmake)

## 2. set up project specific configuration ##
message(STATUS "${Blue}${PROJECT_NAME}-2. Setting up project parameters ...${ColorReset}")

# load custom task
message(STATUS "Load custom task = ${TASK}")
if (NOT ${TASK} STREQUAL "")
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/examples/${TASK}.cmake)
    message(STATUS "Load custom setting - ${TASK} -")
    message(STATUS "Include ${CMAKE_CURRENT_SOURCE_DIR}/cmake/examples/${TASK}.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/examples/${TASK}.cmake)
  else()
    message(FATAL_ERROR "Cannot find custom file cmake/examples/${TASK}.cmake")
  endif()
  set(TASK ${TASK} CACHE STRING "Custom task to be performed")
else()
  set(TASK "default" CACHE STRING "Custom task to be performed")
endif()

message(STATUS "Include ${CMAKE_CURRENT_SOURCE_DIR}/cmake/parameters.cmake")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/parameters.cmake)

## 3. set up project system libraries ##
message(STATUS "${Blue}${PROJECT_NAME}-3. Setting up system libraries ...${ColorReset}")
find_package(Torch REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Disort REQUIRED)
find_package(Harp REQUIRED)
find_package(Kintera REQUIRED)
find_package(Snap REQUIRED)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/fmt.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gtest.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/yamlpp.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/athena.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/application.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/rfm.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/rrtmg.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/pvfmm.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/minichem.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/exofmsrt.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/pybind11.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/c3m.cmake)

## 4. set up project configure file and library ##
message(STATUS "${Blue}${PROJECT_NAME}-4. Setting up project libraries ...${ColorReset}")
configure_file(${CMAKE_SOURCE_DIR}/configure.h.in configure.h @ONLY)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(data)

## 5. set up examples and tests
message(STATUS "${Blue}${PROJECT_NAME}-5. Setting up unit tests ..${ColorReset}")

add_subdirectory(examples)
#add_subdirectory(tests)

## 6. set up python binding ##
#add_subdirectory(python)

## 7. add uninstall targets ##
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(canoe_uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

message(STATUS "${Green}==== ${PROJECT_NAME} configure done ====${ColorReset}")
message(STATUS "")
