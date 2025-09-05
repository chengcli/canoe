# configuration for euler equations on cubed sphere

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

set(MINICHEM ON)
set(NCHEM 13)
# set(NVAPOR 1) set(NCLOUD 2) set(NPHASE_LEGACY 3)
set(NETCDF ON)
set(PNETCDF ON)
set(PT ON)
set(MPI ON)
set(DISORT OFF)
set(PYTHON_BINDINGS ON)
set(RSOLVER hllc)
# set(CUBED_SPHERE ON) set_if_empty(RSOLVER hllc_transform) set(GLOG ON)
