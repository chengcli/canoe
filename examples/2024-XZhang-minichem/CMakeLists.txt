
if (${COORDINATE_SYSTEM} STREQUAL "spherical_polar")
  setup_problem(minichem_sphr)
else()
  setup_problem(minichem)
endif()

# 2. Copy input files to run directory
file(GLOB inputs *.py *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
