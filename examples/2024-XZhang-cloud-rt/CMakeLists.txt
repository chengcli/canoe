
# 1. Compile hywater
#if (${NVAPOR} EQUAL 1 AND ${NCLOUD} EQUAL 2 AND ${NPHASE_LEGACY} EQUAL 3)
if (${NVAPOR} EQUAL 2)
  setup_problem(rtcld)
else()
  # Compile hjupiter if none of the above conditions are met
  setup_problem(rtdry)
endif()

# 2. Copy input files to run directory
file(GLOB inputs *.py *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
