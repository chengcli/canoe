include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  application
  # GIT_REPOSITORY https://github.com/chengcli/application/ GIT_TAG cli/flush)
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/chengcli/application/archive/refs/tags/v0.4.2.tar.gz)

FetchContent_MakeAvailable(application)

include_directories(${application_SOURCE_DIR})
