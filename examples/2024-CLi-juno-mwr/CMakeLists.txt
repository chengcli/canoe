# =====================================================
# Juno Microwave Radiometer Forward and Inversion Model
# =====================================================

# 1. Compile juno_mwr.cpp
setup_problem(juno_mwr)

# 4. Copy input file to run directory
file(GLOB inputs *.inp *.yaml *.py)
foreach(input ${inputs})
  # softlink inp files
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
