
# 1. Compile hywater
if (${NVAPOR} EQUAL 1 AND ${NCLOUD} EQUAL 2 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(hywater)
endif()

# 2. Compile hjupiter
setup_problem(hjupiter)

# 2. Copy input files to run directory
file(GLOB inputs *.py *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
