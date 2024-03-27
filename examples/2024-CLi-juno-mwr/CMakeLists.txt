# =====================================================
# Juno Microwave Radiometer Forward and Inversion Model
# =====================================================

# 1. Compile juno_mwr.cpp
setup_problem(juno_mwr)

# PGEN variable will be used in python bindings
set(PGEN ${CMAKE_CURRENT_SOURCE_DIR}/juno_mwr
  CACHE INTERNAL "Path to the source file for problem generator")

# 4. Copy input file to run directory
file(GLOB inputs *.inp *.yaml *.py)
foreach(input ${inputs})
  # softlink inp files
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
