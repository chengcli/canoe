# ============================
# A collection of unit tests
# ============================

# Enable testing functionality
enable_testing()

if (${HYDROSTATIC})
  add_subdirectory(1d-rad-jupiter)
endif()

# yaml
add_executable(test_yaml_read test_yaml_read.cpp)
target_link_libraries(test_yaml_read gtest_main)
add_test(NAME test_yaml_read COMMAND test_yaml_read)
