# ============================
# A collection of unit tests
# ============================

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# Enable testing functionality
enable_testing()

#if (${HYDROSTATIC})
#  add_subdirectory(1d-rad-jupiter)
#endif()

# yaml
#set(namel test_yaml_read)
#add_executable(${namel}.${buildl} ${namel}.cpp)
#target_link_libraries(${namel}.${buildl} gtest_main yaml-cpp)
#add_test(NAME ${namel}.${buildl} COMMAND ${namel}.${buildl})

file(GLOB src_files *.cpp)

foreach(test ${src_files})
  get_filename_component(name ${test} NAME_WE)

  add_executable(${name}.${buildl} ${name}.cpp)

  set_target_properties(
    ${name}.${buildl}
    PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )

  target_include_directories(${name}.${buildl}
    PRIVATE ${CANOE_INCLUDE_DIR})

  target_link_libraries(${name}.${buildl} gtest_main
    ${CANOE_LIBRARY_${buildu}})

  add_test(NAME ${name}.${buildl} COMMAND ${name}.${buildl})
endforeach()


# copy yaml files
file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()
