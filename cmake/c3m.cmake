# Fetch c3m and build
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set_if_empty(ACCOUNT $ENV{GH_ACCOUNT})
set_if_empty(TOKEN $ENV{GH_TOKEN})

option(C3M "Build C3M" OFF)

if(C3M)
  FetchContent_Declare(
    c3m
    GIT_REPOSITORY https://github.com/chengcli/C3M
    GIT_TAG main)

  FetchContent_MakeAvailable(c3m)
endif()
