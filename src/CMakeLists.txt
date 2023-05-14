# compile the following libraries
#
# climath : this is a math library written in c
# debugger: this is a debugger library written in cpp

add_subdirectory(climath)
set(CLIMATH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"climath include directory")

add_subdirectory(debugger)
set(DEBUGGER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"debugger include directory")
