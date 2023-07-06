# ============================
# A collection of unit tests
# ============================

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

macro(add_test_case test)
  get_filename_component(name ${test} NAME_WE)

  add_executable(${name}.${buildl} ${name}.cpp)

  set_target_properties(
    ${name}.${buildl}
    PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

  target_include_directories(${name}.${buildl}
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CANOE_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${NETCDF_INCLUDES}
    ${PNETCDF_INCLUDE_DIR}
    )

  target_link_libraries(${name}.${buildl} gtest_main
    ${ATHENAPP_LIBRARY_${buildu}}
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
add_test_case(test_absorber.cpp)
add_test_case(test_yaml_read.cpp)

if (${NVAPOR} EQUAL 0)
  add_test_case(test_mesh.cpp)
endif()

if (${NVAPOR} EQUAL 2)
  add_test_case(test_variable.cpp)
  add_test_case(test_radiation.cpp)
  add_test_case(test_microwave_opacity.cpp)
  add_test_case(test_thermodynamics.cpp)
endif()

# copy yaml files
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()
