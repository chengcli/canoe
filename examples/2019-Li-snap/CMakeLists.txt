string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

macro(setup_problem namel)
  set_target_properties(
      ${namel}.${buildl}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
      COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
      )

  target_link_libraries(${namel}.${buildl}
      ${CANOE_LIBRARY_${buildu}}
      ${NETCDF_LIBRARIES}
      ${PNETCDF_LIBRARY}
      ${MPI_CXX_LIBRARIES}
      ${GLOG_LIBRARIES}
      )

  target_include_directories(${namel}.${buildl}
      PRIVATE
      ${CMAKE_BINARY_DIR}
      ${CANOE_INCLUDE_DIR}
      ${EIGEN3_INCLUDE_DIR}
      ${MPI_CXX_INCLUDE_PATH}
      ${NETCDF_INCLUDES}
      ${PNETCDF_INCLUDE_DIR}
      )
endmacro()

# 1. Compile straka and robert
foreach(namel straka robert)
  add_executable(${namel}.${buildl} ${namel}.cpp
    ${CMAKE_BINARY_DIR}/main_${TASK}.cpp)
  setup_problem(${namel})
endforeach()

# 2. Compile bryan
if (${NVAPOR} EQUAL 1 AND ${NCLOUD} EQUAL 1 AND ${NPHASE_LEGACY} EQUAL 2)
  set(namel bryan)

  add_executable(${namel}.${buildl} ${namel}_pgen.cpp ${namel}_mods.cpp
    ${CMAKE_BINARY_DIR}/main_${TASK}.cpp)

  setup_problem(${namel})
endif()

# 3. Compile sedimentation
if (${NVAPOR} EQUAL 2 AND ${NCLOUD} EQUAL 4 AND ${NPHASE_LEGACY} EQUAL 3)
  set(namel sedimentation)

  add_executable(${namel}.${buildl} ${namel}.cpp ../giants_enroll_vapor_functions.cpp
    ${CMAKE_BINARY_DIR}/main_${TASK}.cpp)

  setup_problem(${namel})
endif()

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
