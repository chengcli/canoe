include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  application
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/chengcli/application/archive/refs/tags/v0.7.tar.gz)

FetchContent_MakeAvailable(application)

include_directories(${application_SOURCE_DIR})
