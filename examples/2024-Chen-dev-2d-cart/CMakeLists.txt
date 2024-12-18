# ============================================
# Examples published in Li et al. (2022), PNAS
# ============================================

# 1. Compile shallow_xy and polar_vortex
setup_problem(polar_vortex_cf)

# 1. Copy input file to run directory

file(GLOB inputs *.inp *.dat)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
