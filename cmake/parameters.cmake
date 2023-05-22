## define default parameters

macro(SET_IF_EMPTY _variable)
    if("${${_variable}}" STREQUAL "")
        set(${_variable} ${ARGN})
    endiF()
endmacro()

## populate the default values

SET_IF_EMPTY(NCLOUDS 0)

SET_IF_EMPTY(EQUATION_OF_STATE "ideal_moist")
