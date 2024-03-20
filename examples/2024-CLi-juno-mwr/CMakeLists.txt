# =====================================================
# Juno Microwave Radiometer Forward and Inversion Model
# =====================================================

# 1. Compile juno_mwr.cpp
setup_problem(juno_mwr)

# 4. Copy input file to run directory
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
