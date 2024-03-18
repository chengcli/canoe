
# 1. Compile amars
if (${NVAPOR} EQUAL 4 AND ${NCLOUD} EQUAL 8 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(amars)
endif()

# 2. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
