include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

option(MINICHEM "Build minichem" OFF)

if(MINICHEM)
  FetchContent_Declare(
    mini_chem
    GIT_REPOSITORY https://github.com/chengcli/mini_chem
    GIT_TAG main)
  FetchContent_MakeAvailable(mini_chem)
endif()
