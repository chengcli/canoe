# configuration for euler equations on cubed sphere

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

set(MINICHEM ON)
# set(CUBED_SPHERE ON) set(COORDINATE_SYSTEM "spherical_polar")
set(NVAPOR 1)
set(NCLOUD 2)
set(NPHASE_LEGACY 3)
set(NETCDF ON)
set(PNETCDF ON)
set(MPI ON)
set(DISORT OFF)
set(PYTHON_BINDINGS ON)
set(NCHEMISTRY 12)
set(RSOLVER lmars)
# set_if_empty(RSOLVER hllc_transform) set(GLOG ON)

# set(NUMBER_CHEMICAL_SPECIES 2) set(CHEMISTRY_ENABLED 1) set(CHEMNETWORK_HEADER
# "../chemistry/network/H2.hpp")

# SET_IF_EMPTY(CHEMRADIATION_INTEGRATOR "none")
