# configuration for straka hydrodynamcis

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(PLANET "Earth")
set(CUBED_SPHERE OFF)
set(COORDINATE_SYSTEM "spherical_polar")
set(RSOLVER hllc_transform)
set(TASKLIST ImplicitHydroTasks)
set(NETCDF ON)
set(MPI ON)
set(PNETCDF ON)
