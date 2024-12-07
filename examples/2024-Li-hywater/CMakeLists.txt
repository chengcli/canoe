# =========================
# Examples published in XXX
# =========================

# 1. Compile straka, robert and bryan
setup_problem(hywater)

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
