# configuration for polar vortex simulation

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

set(EOS shallow_xy)
set(PNETCDF ON)
set(MPI ON)
