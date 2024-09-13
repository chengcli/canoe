  setup_problem(minichem)

# 2. Copy input files to run directory
file(GLOB inputs *.py *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
