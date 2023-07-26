# configuration for Mars FarIR forward calculation

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 2)
set_if_empty(NVAPOR 0)

# canoe variables
set_if_empty(NCLOUD 0)
set_if_empty(NTRACER 5)

# canoe task set(TASKLIST InversionTasks)

# canoe configure
set(PLANET "Mars")
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/opt/")
set(HYDROSTATIC ON)
set(RFM ON)
set(NETCDF ON)
set(DISORT ON)
