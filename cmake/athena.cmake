include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

set(decomposition_patch
    git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/19.decomposition.patch)

FetchContent_Declare(
  athenapp
  GIT_REPOSITORY https://github.com/chengcli/athenapp/
  GIT_TAG snap-mods
  PATCH_COMMAND ${decomposition_patch}
  UPDATE_DISCONNECTED TRUE)
# DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL
# https://github.com/chengcli/athenapp/archive/refs/tags/v0.5.tar.gz)

FetchContent_MakeAvailable(athenapp)

include_directories(${athenapp_SOURCE_DIR})
