# =============================================
# Examples published in Li et al. (2018), JQSRT
# =============================================

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# 1. Compile giants_re.cpp

setup_problem(giants_re)

# 2. Copy input file to run directory
file(GLOB inputs *.inp *.toml *.yaml)
foreach(input ${inputs})
  file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
