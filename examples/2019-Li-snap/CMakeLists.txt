# ============================================
# Examples published in Li et al. (2019), ApJS
# ============================================

# 1. Compile straka, robert and bryan
setup_problem(straka)
setup_problem(robert)
setup_problem(bryan)

# 2. Compile sedimentation
if (${NVAPOR} EQUAL 2 AND ${NCLOUD} EQUAL 4 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(sedimentation)
  setup_problem(jvortex)
endif()

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
