# configure file for test sedimentation

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(NVAPOR 2)
set(NCLOUD 4)
set(NPHASE_LEGACY 3)
set(NETCDF ON)
set(PNETCDF ON)
set(MPI ON)
set(TASKLIST ImplicitHydroTasks)
set(RSOLVER lmars)
