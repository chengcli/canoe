# configuration for Jupiter Juno forward and inversion calculation

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 2)
set_if_empty(NVAPOR 0)

# canoe variables
set_if_empty(NCLOUD 0)
set_if_empty(NTRACER 0)

# canoe task set(TASKLIST InversionTasks)

# canoe configure
set(PLANET "Uranus")
set(HYDROSTATIC ON)
set(RFM ON)
set(NETCDF ON)
set(DISORT ON)
