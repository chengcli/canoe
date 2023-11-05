# ============================================
# Examples published in Li et al. XXX
# ============================================

# 1. Compile uranus_mwr.cpp
setup_problem(uranus_mwr)

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
