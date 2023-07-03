# configuration for robert hydrodynamcis

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(PLANET "Earth")
set(MPI ON)
set(PNETCDF ON)
set(TASKLIST ImplicitHydroTasks)
