# ===============================================
# Examples published in Cheng et al. (2023), ApJS
# ===============================================

# 1. Compile dam_break.cpp
setup_problem(dam_break)

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.dat)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
