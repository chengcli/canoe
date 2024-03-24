include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

option(RRTMG_SW "Build RRTMG_SW" OFF)

if(RRTMG_SW)
  FetchContent_Declare(
    RRTMG_SW
    GIT_REPOSITORY https://github.com/chengcli/RRTMG_SW
    GIT_TAG master)
  FetchContent_MakeAvailable(RRTMG_SW)
endif()

option(RRTMG_LW "Build RRTMG_LW" OFF)

if(RRTMG_LW)
  FetchContent_Declare(
    RRTMG_LW
    GIT_REPOSITORY https://github.com/chengcli/RRTMG_LW
    GIT_TAG master)
  FetchContent_MakeAvailable(RRTMG_LW)
endif()
