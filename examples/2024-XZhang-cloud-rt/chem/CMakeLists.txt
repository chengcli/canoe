
# 1. Compile hywater
if (${NVAPOR} EQUAL 1 AND ${NCLOUD} EQUAL 2 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(hywater)
endif()


set(namel hjupiter)
add_executable(${namel}.${buildl} ${namel}.cpp H2.cpp ${CMAKE_SOURCE_DIR}/main.cpp)

set_target_properties(
${namel}.${buildl}
PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
			 COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}})

target_link_libraries(${namel}.${buildl} ${CANOE_LIBRARY_${buildu}})

target_include_directories(
${namel}.${buildl}
PRIVATE ${CMAKE_BINARY_DIR} ${CANOE_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}
		${MPI_CXX_INCLUDE_PATH} ${NETCDF_INCLUDES} ${PNETCDF_INCLUDE_DIR})

# 2. Copy input files to run directory
file(GLOB inputs *.py *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
