# configuration for euler equations on cubed sphere

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

set(CUBED_SPHERE ON)
set(NETCDF OFF)
set(PNETCDF ON)
set(MPI ON)
set(DISORT ON)
set(PYTHON_BINDINGS ON)
set(NCHEMISTRY 1)
set_if_empty(RSOLVER hllc_transform)
# set(GLOG ON)

#set(NUMBER_CHEMICAL_SPECIES 2)
#set(CHEMISTRY_ENABLED 1)
#set(CHEMNETWORK_HEADER "../chemistry/network/H2.hpp")

#SET_IF_EMPTY(CHEMRADIATION_INTEGRATOR "none")
