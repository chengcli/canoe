set(namel c3m)
string(TOUPPER ${namel} nameu)

file(GLOB src_files
    *.cpp
    )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
    OBJECT
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
    )


target_include_directories(${namel}_${buildl}
    SYSTEM PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${CANTERA_INCLUDE_DIR}
    ${ATHENAPP_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}
    )

target_link_libraries(${namel}_${buildl}
    PRIVATE
    ${CANTERA_LIBRARY}
    ${ATHENAPP_LIBRARY_${build}}
    python3.9
    dl
    )


