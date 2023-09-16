# configuration for euler equations on cubed sphere

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(CUBED_SPHERE ON)
set_if_empty(RSOLVER hllc_transform)
set(TASKLIST ImplicitHydroTasks)
set(NETCDF ON)
set(MPI ON)
set(PNETCDF ON)
