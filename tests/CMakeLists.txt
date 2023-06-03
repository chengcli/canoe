# ============================
# A collection of unit tests
# ============================

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)

# Enable testing functionality
enable_testing()

#if (${HYDROSTATIC})
#  add_subdirectory(1d-rad-jupiter)
#endif()

# yaml
set(namel test_yaml_read)
add_executable(${namel}.${buildl} ${namel}.cpp)
target_link_libraries(${namel}.${buildl} gtest_main yaml-cpp)
add_test(NAME ${namel}.${buildl} COMMAND ${namel}.${buildl})

# copy yaml files
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()
