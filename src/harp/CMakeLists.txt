# Installs harp library
#
#   library: libharp_debug.o
#   library: libharp_release.o
#
# Define the following variables
#
#   HARP_INCLUDE_DIR
#   HARP_LIBRARY_DEBUG
#   HARP_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${INVERSION_INCLUDE_DIR})
#   target_link_libraries( ${INVERSION_LIBRARY_DEBUG})

set(namel harp)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
    *.cpp
    )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
    OBJECT
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

if(DISORT)
  target_include_directories(${namel}_${buildl}
      PRIVATE
      ${NETCDF_INCLUDES}
      ${TOMLPP_INCLUDE_DIR}
      ${CDISORT_INCLUDE_DIR}
      ${PYTHON_INCLUDE_DIR}
      ${PYBIND11_INCLUDE_DIR}
      ${PYDISORT_INCLUDE_DIR}
      ${CANTERA_INCLUDE_DIR}/cantera/ext  # for yaml-cpp
      )
else()
  target_include_directories(${namel}_${buildl}
      PRIVATE
      ${NETCDF_INCLUDES}
      ${CANTERA_INCLUDE_DIR}/cantera/ext  # for yaml-cpp
      )
endif()
