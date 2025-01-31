# ============================================
# Examples published in Li et al. (2019), ApJS
# ============================================

# 1. Compile straka, robert and bryan
setup_problem(titan)

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
