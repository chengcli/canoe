# configuration for ktable test

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 0)
set_if_empty(NVAPOR 1)

# canoe variables
set_if_empty(NCLOUD 0)
set_if_empty(NTRACER 2)

# canoe configure
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/opt/")
# set(HYDROSTATIC OFF)
set(RFM ON)
set(NETCDF ON)
set(PYTHON_BINDINGS ON)
set(DISORT ON)
