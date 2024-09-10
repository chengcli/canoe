# ============================================
# Examples published in Zhang et al. XXX
# ============================================

# 1. Compile expansion.cpp
setup_problem(expansion)

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
