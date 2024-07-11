
if (${NVAPOR} EQUAL 2 AND ${NCLOUD} EQUAL 4 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(jupiter_crm)
elseif (${NVAPOR} EQUAL 0)
  setup_problem(jupiter_dry)
elseif (${NVAPOR} EQUAL 1)
  setup_problem(jupiter_h2o)
endif()

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
