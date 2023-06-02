include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  athenapp
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/chengcli/athenapp/archive/refs/tags/v0.3.1-alpha.tar.gz
)

FetchContent_GetProperties(athenapp)

if(NOT athenapp_POPULATED)
  FetchContent_Populate(athenapp)
  add_subdirectory(${athenapp_SOURCE_DIR} ${athenapp_BINARY_DIR})
endif()
