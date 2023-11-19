# configure file for hydroge-water world

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(NVAPOR 1)
set(NCLOUD 1)
set(NPHASE_LEGACY 2)
set(NETCDF ON)
set(PNETCDF ON)
set(MPI ON)
set(DISORT ON)
set(PYTHON_BINDINGS OFF)
set(TASKLIST ImplicitHydroTasks)
set(RSOLVER lmars)
