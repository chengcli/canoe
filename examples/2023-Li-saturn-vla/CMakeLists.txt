# ========================================================
# Examples published in Li et al. (2019), Science Advances
# ========================================================

# 1. Compile saturn_radio.cpp
setup_problem(saturn_radio)

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
