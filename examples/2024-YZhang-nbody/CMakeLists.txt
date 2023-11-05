# ============================================
# Examples published in Zhang et al. XXX
# ============================================

# 1. Compile nbody.cpp
setup_problem(nbody)

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
