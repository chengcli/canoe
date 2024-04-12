# configure file for test jupiter crm

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(NVAPOR 1)
set(NCLOUD 2)
set(NPHASE_LEGACY 3)
set(PNETCDF ON)
set(MPI ON)
set(TASKLIST ImplicitHydroTasks)
set(RSOLVER lmars)
