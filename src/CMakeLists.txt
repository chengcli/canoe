# compile the following modules
#
# climath : math library in c
# debugger: debugger library in cpp
# astro: astronomy library in cpp
# utils: utility library in c/cpp
# harp: radiation library in c/cpp (backed by cdisort)
# inversion: spectral inversion library
# outputs: write netcdf/pnetcdf outputs

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
  )

string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

#set(CANOE_LIBRARY_${buildu} "" CACHE "CANOE LIBRARIES")
set(CANOE_LIBRARY_${buildu} "item1" "item2" "item3" CACHE INTERNAL "A list")

add_subdirectory(climath)
add_subdirectory(debugger)
add_subdirectory(astro)
add_subdirectory(utils)
add_subdirectory(snap)
add_subdirectory(harp)
add_subdirectory(inversion)
add_subdirectory(outputs)
