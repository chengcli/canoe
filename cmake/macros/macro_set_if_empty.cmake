# A small macro used for set a variable if it is empty
#
# Usage: set_if_empty(name)

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()
