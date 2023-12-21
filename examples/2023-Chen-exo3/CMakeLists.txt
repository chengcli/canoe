# ===============================================
# Examples published in Cheng et al. (2023), ApJS
# ===============================================

# 1. Compile dam_break.cpp
if(${HYDROSTATIC})
  if(${EQUATION_OF_STATE} STREQUAL "shallow_yz")
    setup_problem(dam_break)
    setup_problem(test_blast)
    setup_problem(W92)
    setup_problem(test_injection)
    setup_problem(steady_zonal)
  endif()
else()
  setup_problem(test_adiabat)
  setup_problem(hs94)
  setup_problem(hot_jupiter)
  setup_problem(polar_dry)
endif()

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.dat)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
