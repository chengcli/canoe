# ============================
# A collection of unit tests
# ============================

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

configure_file(globals.cpp.in globals.cpp @ONLY)

macro(add_test_giants test)
  get_filename_component(name ${test} NAME_WE)

  add_executable(${name}.${buildl} ${name}.cpp
	  globals.cpp
    ${CMAKE_SOURCE_DIR}/examples/giants_add_absorber.cpp
    ${CMAKE_SOURCE_DIR}/examples/giants_enroll_vapor_functions.cpp
	  )

  set_target_properties(
    ${name}.${buildl}
    PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

  target_include_directories(${name}.${buildl}
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CANOE_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${NETCDF_INCLUDES}
    ${PNETCDF_INCLUDE_DIR}
    )

  target_link_libraries(${name}.${buildl} gtest_main
    ${CANOE_LIBRARY_${buildu}}
    ${NETCDF_LIBRARIES}
    ${PNETCDF_LIBRARY}
    ${MPI_CXX_LIBRARIES}
    )

  add_test(NAME ${name}.${buildl} COMMAND ${name}.${buildl})
endmacro()

macro(add_test_case test)
  get_filename_component(name ${test} NAME_WE)

  add_executable(${name}.${buildl} ${name}.cpp globals.cpp)

  set_target_properties(
    ${name}.${buildl}
    PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

  target_include_directories(${name}.${buildl}
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CANOE_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${NETCDF_INCLUDES}
    ${PNETCDF_INCLUDE_DIR}
    )

  target_link_libraries(${name}.${buildl} gtest_main
    ${CANOE_LIBRARY_${buildu}}
    ${NETCDF_LIBRARIES}
    ${PNETCDF_LIBRARY}
    ${MPI_CXX_LIBRARIES}
    )

  add_test(NAME ${name}.${buildl} COMMAND ${name}.${buildl})
endmacro()

# Enable testing functionality
enable_testing()

# Add test cases
add_test_case(test_weno.cpp)
add_test_giants(test_absorber.cpp)
add_test_giants(test_yaml_read.cpp)

if (${NVAPOR} EQUAL 0)
  add_test_case(test_mesh.cpp)
endif()

if (${NVAPOR} EQUAL 1 AND ${NCLOUD} EQUAL 1 AND ${NPHASE_LEGACY} EQUAL 2)
  add_test_case(test_moist_adiabat.cpp)
endif()

if (${NVAPOR} EQUAL 2)
  if (${NCLOUD} EQUAL 4)
    add_test_giants(test_impl.cpp)
    add_test_giants(test_variable.cpp)
    add_test_giants(test_radiation.cpp)
    add_test_giants(test_microwave_opacity.cpp)
    add_test_giants(test_thermodynamics.cpp)
  elseif (${NCLOUD} EQUAL 5)
    add_test_giants(test_ammonium_hydrosulfide.cpp)
  endif()
endif()

# copy yaml files
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()
