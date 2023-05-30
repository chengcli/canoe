string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# 1. Compile test_jupiter_re.cpp

set(namel test_jupiter_re)

add_executable(${namel}.${buildl} ${namel}.cpp)

set_target_properties(
    ${namel}.${buildl}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

target_link_libraries(test_jupiter_re.${buildl}
    ${CPPDISORT_LIBRARY_${buildu}}
    ${ATHENA_LIBRARY_${buildu}}
    ${HARP_LIBRARY_${buildu}}
    )

target_include_directories(
  ${namel}.${buildl}
  PRIVATE
  ${ATHENAPP_INCLUDE_DIR}
  ${MODULES_INCLUDE_DIR}
  ${CPPDISORT_INCLUDE_DIR}
  ${CMAKE_BINARY_DIR}
  ${TOMLPP_INCLUDE_DIR}
  )

# 1. Copy input file to run directory

file(GLOB inputs *.inp *.toml)
foreach(input ${inputs})
    file(COPY input DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
