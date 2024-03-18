# ==========================
# Examples published in XXXX
# ==========================

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# 1. Compile jupiter1d.cpp
setup_problem(jupiter1d)

# 2. Copy input file to run directory
file(GLOB inputs *.py *.yaml *.inp)
foreach(input ${inputs})
  # softlink inp files
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
