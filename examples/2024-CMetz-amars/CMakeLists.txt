
# 1. Compile amars
if (${NVAPOR} EQUAL 3 AND ${NCLOUD} EQUAL 6 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(amars)
endif()

if (${NVAPOR} EQUAL 1 AND ${NCLOUD} EQUAL 2 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(amarsw)
endif()

# 2. Copy input files to run directory
file(GLOB inputs *.inp *.yaml *.py)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
