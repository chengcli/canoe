## set up transport tests ##

file(GLOB _tests *.cc *cpp)
foreach(_test ${_tests})
  foreach(build ${BUILD_TYPES})
    string(TOLOWER ${build} buildl)
    get_filename_component(_name ${_test} NAME_WE)

    add_executable(${_name}.${buildl} ${_test})

    set_target_properties(${_name}.${buildl}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${EXECUTABLE_RELDIR}"
      COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${build}}
      )

    message(STATUS ${ICECHEM_INCLUDE_DIUR})

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      # Intel (at least up to 19.0.5) confuses an external TBB installation with
      # our selected one if we use -Isystem includes. Work around this by using
      # normal includes.
      # See https://github.com/dealii/dealii/issues/8374 for details.
      target_include_directories(${_name}.${buildl}
        PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${DEAL_II_INCLUDE_DIRS}
        ${ATHENAPP_INCLUDE_DIR}
        ${ICECHEM_INCLUDE_DIR}
        ${CMAKE_BINARY_DIR}
        )
    else()
      target_include_directories(${_name}.${buildl}
        SYSTEM PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${DEAL_II_INCLUDE_DIRS}
        ${ATHENAPP_INCLUDE_DIR}
        ${ICECHEM_INCLUDE_DIR}
        ${CMAKE_BINARY_DIR}
        )
    endif()

    target_link_libraries(${_name}.${buildl}
      ${DEAL_II_TARGET_RELEASE}
      ${ICECHEM_LIBRARY_${build}}
      )
  endforeach()
endforeach()

## copy plot subroutines to executable folder ##
file(GLOB _inputs *.py)
foreach(_input ${_inputs})
  get_filename_component(_name ${_input} NAME)

  configure_file(${_input} ${CMAKE_BINARY_DIR}/${EXECUTABLE_RELDIR}/${_name} COPYONLY)
endforeach()
