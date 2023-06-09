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

target_link_libraries(${namel}.${buildl}
    ${ATHENAPP_LIBRARY_${buildu}}
    ${CANOE_LIBRARY_${buildu}}
    ${NETCDF_LIBRARIES}
    ${CPPDISORT_LIBRARY_${buildu}}
    )

target_include_directories(${namel}.${buildl}
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src
    )

# 1. Copy input file to run directory

file(GLOB inputs *.inp *.toml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
