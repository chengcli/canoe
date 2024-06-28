# configure file for rad unit test case

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
set(NTRACER 2)
# set(NPHASE_LEGACY 2)
set(NETCDF ON)
set(PNETCDF ON)
set(MPI ON)
set(RSOLVER lmars)
set(DISORT ON)
set(PYTHON_BINDINGS ON)
set(RFM ON)
