# ==========================
# Examples published in XXXX
# ==========================

# 2. Copy input file to run directory
file(GLOB inputs *.py *.yaml)
foreach(input ${inputs})
  # softlink inp files
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
