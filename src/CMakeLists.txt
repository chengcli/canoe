# compile the following modules
#
# climath : math library in c
# astro: astronomy library in cpp
# utils: utility library in c/cpp
# outputs: write netcdf/pnetcdf outputs

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
  ${MPI_CXX_INCLUDE_PATH}
  ${MPI_CXX_HEADER_DIR}
  ${NETCDF_INCLUDES}
  ${PNETCDF_INCLUDE_DIR}
  SYSTEM ${TORCH_API_INCLUDE_DIR}
  SYSTEM ${TORCH_INCLUDE_DIR}
  SYSTEM ${EIGEN3_INCLUDE_DIR}
  SYSTEM ${DISORT_INCLUDE_DIR}
  SYSTEM ${HARP_INCLUDE_DIR}
  SYSTEM ${KINTERA_INCLUDE_DIR}
  SYSTEM ${SNAP_INCLUDE_DIR}
  SYSTEM ${FMT_INCLUDE_DIR}
  )

set(namel canoe)
string(TOUPPER ${namel} nameu)

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

file(GLOB src_files
    *.cpp
    climath/*.c
    astro/*.cpp
    astro/Jupiter/*.cpp
    astro/Saturn/*.cpp
    utils/*.cpp
    tasklist/*.cpp
    exo3/*.cpp
    exchanger/*.cpp
    )

add_library(${namel}_${buildl}
    OBJECT
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

target_link_libraries(${namel}_${buildl}
    PUBLIC
    ${NETCDF_LIBRARIES}
    ${PNETCDF_LIBRARY}
    ${TORCH_LIBRARY}
    ${TORCH_CPU_LIBRARY}
    ${C10_LIBRARY}
    ${DISORT_LIBRARY}
    ${HARP_LIBRARY}
    ${KINTERA_LIBRARY}
    ${SNAP_LIBRARY}
    ${MPI_CXX_LIBRARIES}
    ${GLOG_LIBRARIES}
    ${OpenMP_libomp_LIBRARY}
    ${EXTRA_LIBS}
    fmt::fmt
    yaml-cpp::yaml-cpp
    )

add_subdirectory(snap)
add_subdirectory(outputs)

if(UNIX AND NOT APPLE)
  set(EXTRA_LIBS dl stdc++fs)
endif()

set(CANOE_LIBRARY_${buildu}
  "canoe_${buildl}"
  "athenapp_${buildl}"
  "application_${buildl}"
  "snap_${buildl}"
  "outputs_${buildl}"
  $<$<BOOL:${MINICHEM}>:minichem::minichem++>
  $<$<BOOL:${EXOFMSRT}>:exofmsrt::exofmsrt++>
  CACHE INTERNAL "canoe library")

set(CANOE_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}"
  CACHE INTERNAL "canoe include directory")
