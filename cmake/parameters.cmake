# define default parameters

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# populate the default values

set_if_empty(NCLOUDS 0)

set_if_empty(EQUATION_OF_STATE "ideal_moist")
