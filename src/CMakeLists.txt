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
  )

set(namel canoe)
string(TOUPPER ${namel} nameu)

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

file(GLOB src_files *.cpp)

add_library(${namel}_${buildl}
    OBJECT
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

add_subdirectory(climath)
add_subdirectory(astro)
add_subdirectory(utils)
add_subdirectory(snap)
add_subdirectory(outputs)
add_subdirectory(tasklist)
add_subdirectory(exo3)
add_subdirectory(exchanger)
add_subdirectory(special)

if(UNIX AND NOT APPLE)
  set(EXTRA_LIBS dl stdc++fs)
endif()

set(CANOE_LIBRARY_${buildu}
  "canoe_${buildl}"
  "athenapp_${buildl}"
  "application_${buildl}"
  "climath_${buildl}"
  "astro_${buildl}"
  "utils_${buildl}"
  "snap_${buildl}"
  "outputs_${buildl}"
  "dusts_${buildl}"
  "flask_${buildl}"
  "tracer_${buildl}"
  "tasklist_${buildl}"
  "exo3_${buildl}"
  "exchanger_${buildl}"
  $<$<BOOL:${MINICHEM}>:minichem::minichem++>
  $<$<BOOL:${EXOFMSRT}>:exofmsrt::exofmsrt++>
  ${YAML_CPP_LIBRARIES}
  ${CPPDISORT_LIBRARY_${buildu}}
  ${PYTHON_LIBRARY_RELEASE}
  ${NETCDF_LIBRARIES}
  ${PNETCDF_LIBRARY}
  ${MPI_CXX_LIBRARIES}
  ${GLOG_LIBRARIES}
  ${OpenMP_libomp_LIBRARY}
  ${CANTERA_LIBRARY}
  ${TORCH_LIBRARY}
  ${TORCH_CPU_LIBRARY}
  ${C10_LIBRARY}
  ${DISORT_LIBRARY}
  ${HARP_LIBRARY}
  ${KINTERA_LIBRARY}
  ${EXTRA_LIBS}
  CACHE INTERNAL "canoe library")

set(CANOE_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}"
  CACHE INTERNAL "canoe include directory")
