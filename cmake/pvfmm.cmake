# Fetch pvfmm and build
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

option(PVFMM "Build PVFMM" OFF)

if (PVFMM)
  set(patch_command
      git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/pvfmm.01.sctl_stracktrace.patch
      )

  FetchContent_Declare(
    PVFMM
    GIT_REPOSITORY https://github.com/chengcli/pvfmm/
    GIT_TAG master
    PATCH_COMMAND ${patch_command}
    UPDATE_DISCONNECTED TRUE)

  FetchContent_MakeAvailable(PVFMM)

  include_directories(${PVFMM_SOURCE_DIR} ${PVFMM_BINARY_DIR})
endif()
