# configure file for hydrogen-water world

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

set(CUBED_SPHERE ON) 
#set(NBLOCKS 600)
#set(COORDINATE_SYSTEM "spherical_polar")
set(NVAPOR 1)
set(NCLOUD 2)
set(NPHASE_LEGACY 3)
set(NETCDF ON)
#set(PNETCDF OFF)
set(MPI ON)
set(DISORT ON)
set(PYTHON_BINDINGS ON)
# set_if_empty(RSOLVER hllc_transform) 
# set(GLOG ON)
set(TASKLIST ImplicitHydroTasks) 
set(RSOLVER hllc_transform)
#set(RSOLVER lmars)
