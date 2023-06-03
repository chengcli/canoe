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
message(STATUS "Include ${CMAKE_CURRENT_SOURCE_DIR}/cmake/parameters.cmake")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/parameters.cmake)

## 3. set up project system libraries ##
message(STATUS "3. Set up system libraries")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/athena.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/eigen.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/yamlpp.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/athena.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/rfm.cmake)

## 4. set up project configure file and library ##
message(STATUS "4. Set up project libraries")
configure_file(${CMAKE_SOURCE_DIR}/configure.hpp.in configure.hpp @ONLY)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(data)

## 5. set up python binding ##
message(STATUS "5. Set up python binding")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/pybind11.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/pybind11.cmake)
#add_subdirectory(python)

## 6. set up unit tests
message(STATUS "6. Set up unit tests")
add_subdirectory(tests)

## 7. install library and executables ##
message(STATUS "7. Install python library and executables")
