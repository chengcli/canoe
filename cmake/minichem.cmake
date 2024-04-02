include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

option(MINICHEM "Build minichem" OFF)

if(MINICHEM)
  FetchContent_Declare(
    minichem
    GIT_REPOSITORY https://github.com/chengcli/mini_chem
    GIT_TAG main)
  FetchContent_MakeAvailable(minichem)

  include_directories(${minichem_SOURCE_DIR})
endif()

