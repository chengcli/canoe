# ==========================
# Examples published in XXXX
# ==========================

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# 2. Copy input file to run directory
file(GLOB inputs *.py *.yaml *.inp)
foreach(input ${inputs})
  # softlink inp files
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
