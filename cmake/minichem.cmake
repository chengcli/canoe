include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

option(MINICHEM "Build minichem" OFF)

if(MINICHEM)
set(patch_command
    git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/35.mini_ch_i_dlsode.f90.patch)

  FetchContent_Declare(
    minichem
    GIT_REPOSITORY https://github.com/chengcli/mini_chem
    GIT_TAG main
    PATCH_COMMAND ${patch_command}
    UPDATE_DISCONNECTED TRUE)
  FetchContent_MakeAvailable(minichem)

  include_directories(${minichem_SOURCE_DIR})
endif()

