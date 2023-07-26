# Installs utility library
#
#   python: combine.py
#   python: main2mcmc.py
#   python: write_workflow.py
#   c++: fitsmerge.release

message(STATUS "set up utilities ...")

add_executable(kcoeff.release kcoeff.cpp)

set_target_properties(kcoeff.release
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELEASE}
  )

target_include_directories(kcoeff.release
  SYSTEM PRIVATE
  ${CMAKE_BINARY_DIR}
  ${NETCDF_INCLUDES}
  )

target_link_libraries(kcoeff.release
  ${NETCDF_LIBRARIES}
  )

add_executable(hitbin.release hitbin.f90)

set_target_properties(hitbin.release
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE}
  )

file(GLOB pyfiles *.py)

foreach(pyfile ${pyfiles})
  file(COPY ${pyfile}
    DESTINATION
    ${CMAKE_BINARY_DIR}/bin
  )
endforeach()

file(COPY clean_rfm_runs.sh DESTINATION ${CMAKE_BINARY_DIR}/bin)

configure_file(run_ktable.py.in ${CMAKE_BINARY_DIR}/bin/run_ktable.py @ONLY)
configure_file(cktable.py.in ${CMAKE_BINARY_DIR}/bin/cktable.py @ONLY)

# copy yaml and inp files
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin/)
endforeach()
