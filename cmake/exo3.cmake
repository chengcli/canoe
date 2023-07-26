# configuration for straka hydrodynamcis

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 2)

# canoe configure
set(PLANET "Earth")
set(EOS "shallow_yz")
set(HYDROSTATIC ON)
set(CUBED_SPHERE ON)
set(NETCDF ON)
set(MPI ON)
set(PNETCDF ON)
