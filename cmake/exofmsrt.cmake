include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

option(EXOFMSRT "Build exofms-ck" OFF)

if(EXOFMSRT)
  FetchContent_Declare(
    exofms-ck
    GIT_REPOSITORY https://github.com/chengcli/Exo-FMS_column_ck
    GIT_TAG main)
  FetchContent_MakeAvailable(exofms-ck)

  include_directories(${exofmsck_SOURCE_DIR})
endif()
