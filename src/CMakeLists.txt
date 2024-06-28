# compile the following modules
#
# climath : math library in c
# astro: astronomy library in cpp
# utils: utility library in c/cpp
# harp: radiation library in c/cpp (backed by cdisort)
# inversion: spectral inversion library
# outputs: write netcdf/pnetcdf outputs
# opacity: opacity library

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
  ${NETCDF_INCLUDES}
  ${CFITSIO_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_PATH}
  ${MPI_CXX_HEADER_DIR}
  ${PNETCDF_INCLUDE_DIR}
  ${CANTERA_INCLUDE_DIR}
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

target_include_directories(${namel}_${buildl}
    PRIVATE
    ${CANTERA_INCLUDE_DIR}/cantera/ext  # for yaml-cpp
    )

add_subdirectory(climath)
add_subdirectory(astro)
add_subdirectory(utils)
add_subdirectory(snap)
#add_subdirectory(harp)
#add_subdirectory(inversion)
add_subdirectory(outputs)
#add_subdirectory(opacity)
#add_subdirectory(microphysics)
#add_subdirectory(flask)
#add_subdirectory(tracer)
add_subdirectory(tasklist)
add_subdirectory(exo3)
add_subdirectory(exchanger)
#add_subdirectory(special)
#add_subdirectory(single_column)
add_subdirectory(diagnostics)
add_subdirectory(forcing)
#add_subdirectory(nbody)

set(CANOE_LIBRARY_${buildu}
  "canoe_${buildl}"
  "athenapp_${buildl}"
  "application_${buildl}"
  "climath_${buildl}"
  "astro_${buildl}"
  "utils_${buildl}"
  "snap_${buildl}"
  "harp_${buildl}"
  "inversion_${buildl}"
  "outputs_${buildl}"
  "opacity_${buildl}"
  "dusts_${buildl}"
  "flask_${buildl}"
  "tracer_${buildl}"
  "tasklist_${buildl}"
  "exo3_${buildl}"
  "exchanger_${buildl}"
  #"nbody_${buildl}"
  "scm_${buildl}"
  "diagnostics_${buildl}"
  "forcing_${buildl}"
  #"yaml-cpp"
  ${CPPDISORT_LIBRARY_${buildu}}
  ${PYTHON_LIBRARY_RELEASE}
  ${NETCDF_LIBRARIES}
  ${PNETCDF_LIBRARY}
  ${CFITSIO_LIBRARY}
  ${MPI_CXX_LIBRARIES}
  ${GLOG_LIBRARIES}
  ${OpenMP_libomp_LIBRARY}
  ${CANTERA_LIBRARY}
  CACHE INTERNAL "canoe library")

set(CANOE_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}"
  CACHE INTERNAL "canoe include directory")
