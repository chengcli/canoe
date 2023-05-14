# compile the following libraries
#
# climath : this is a math library written in C

add_subdirectory(climath)
set(CLIMATH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	CACHE PATH 
	"climath include directory")
