# compile the following libraries
#
# climath : math library in c
# debugger: debugger library in cpp
# astro: astronomy library in cpp
# utils: utility library in c/cpp
# harp: radiation library in c/cpp (backend by cdisort)

add_subdirectory(climath)
set(CLIMATH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"climath include directory")

add_subdirectory(debugger)
set(DEBUGGER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"debugger include directory")

add_subdirectory(astro)
set(ASTRO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"astro include directory")

add_subdirectory(utils)
set(UTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"utils include directory")
