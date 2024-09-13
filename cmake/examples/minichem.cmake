# configuration for euler equations on cubed sphere

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

set(MINICHEM ON)
set(NCHEMISTRY 13)
set(NVAPOR 1)
set(NCLOUD 2)
set(NPHASE_LEGACY 3)
set(NETCDF OFF)
set(PNETCDF ON)
set(MPI ON)
set(DISORT OFF)
set(PYTHON_BINDINGS ON)
set(RSOLVER lmars)
# set(CUBED_SPHERE ON) set_if_empty(RSOLVER hllc_transform) set(GLOG ON)
