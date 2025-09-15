set(namel outputs)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
  combine_blocks.cpp
  mppnccombine.cpp
  load_user_output_data.cpp
  output_utils.cpp
  debug.cpp
  netcdf.cpp
  pnetcdf.cpp
  fits.cpp
  pt.cpp
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

target_include_directories(${namel}_${buildl}
    PRIVATE
    ${NETCDF_INCLUDES}
    )
