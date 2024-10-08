# =============================================
# Examples published in Li et al. (2018), JQSRT
# =============================================

# 1. Compile giants_re.cpp

setup_problem(giants_re)

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.yaml *.py *.atm)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
