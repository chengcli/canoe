# A small macro used for setting up the build of a test.
#
# Usage: setup_test(name)

macro(setup_test namel)
  string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
  string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

  add_executable(${namel}.${buildl} ${namel}.cpp)

  set_target_properties(${namel}.${buildl}
                        PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}})

  target_include_directories(
    ${namel}.${buildl}
    PRIVATE ${CMAKE_BINARY_DIR}
            ${CANOE_INCLUDE_DIR}
            ${EIGEN3_INCLUDE_DIR}
            ${MPI_CXX_INCLUDE_PATH}
            ${MPI_CXX_HEADER_DIR}
            ${NETCDF_INCLUDES}
            ${PNETCDF_INCLUDE_DIR}
            ${OpenMP_CXX_INCLUDE_DIR}
            ${FFTW_INCLUDE_DIRS}
            ${CANTERA_INCLUDE_DIR}
            ${TORCH_INCLUDE_DIR}
            ${TORCH_API_INCLUDE_DIR})

  target_link_libraries(${namel}.${buildl} gtest_main
                        ${CANOE_LIBRARY_${buildu}})

  add_test(NAME ${namel}.${buildl} COMMAND ${namel}.${buildl})
endmacro()
