# 1. Compile
setup_problem(geyser_wall)

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml *.cpp *.hpp)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
