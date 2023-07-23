# configuration for Jupiter Juno forward and inversion calculation

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 2)
set_if_empty(NVAPOR 4)

# canoe variables
set_if_empty(NCLOUD 9)
set_if_empty(NTRACER 2)

# canoe task set(TASKLIST InversionTasks)

# canoe configure
set(PLANET "Uranus")
set(HYDROSTATIC ON)
set(NETCDF ON)
set(FITS ON)
set(RadianceSolver lambert)
