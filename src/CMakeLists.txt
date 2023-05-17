# compile the following modules
#
# climath : math library in c
# debugger: debugger library in cpp
# astro: astronomy library in cpp
# utils: utility library in c/cpp
# harp: radiation library in c/cpp (backed by cdisort)
# inversion: spectral inversion library
# outputs: write netcdf/pnetcdf outputs

set(MODULES_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"canoe module include directory")

add_subdirectory(astro)
add_subdirectory(climath)
add_subdirectory(debugger)
add_subdirectory(utils)
add_subdirectory(harp)
add_subdirectory(inversion)
add_subdirectory(outputs)
