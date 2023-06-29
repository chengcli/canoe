# configuration for Jupiter Juno forward and inversion calculation

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 2)
set_if_empty(NVAPOR_VARIABLES 3)

# canoe variables
set_if_empty(NCLOUDS 5)
set_if_empty(NTRACER 3)
set_if_empty(NCHEMISTRY 0)

# canoe configure
set(HYDROSTATIC ON)
set(NETCDF ON)
set(FITS ON)
set(RadianceSolver lambert)
