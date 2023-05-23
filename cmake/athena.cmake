include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  athenapp
  URL https://github.com/chengcli/athenapp/archive/refs/tags/v0.1.3-alpha.tar.gz
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE)

FetchContent_GetProperties(athenapp)

if(NOT athenapp_POPULATED)
  FetchContent_Populate(athenapp)
  add_subdirectory(${athenapp_SOURCE_DIR} ${athenapp_BINARY_DIR})
endif()
