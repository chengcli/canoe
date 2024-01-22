# ============================================
# Examples published in Li et al. (2019), ApJS
# ============================================

# 1. Compile straka, robert and bryan
setup_problem(straka)
setup_problem(lstraka)
setup_problem(robert)
setup_problem(bryan)

# 2. Compile sedimentation
if (${NVAPOR} EQUAL 2 AND ${NCLOUD} EQUAL 4 AND ${NPHASE_LEGACY} EQUAL 3)
  setup_problem(sedimentation)
endif()

# 3. Compile jupiter crm
setup_problem(jupiter_crm)

# 4. Copy input files to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
