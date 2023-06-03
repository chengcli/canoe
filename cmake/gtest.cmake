include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

FetchContent_Declare(
  gtest
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.tar.gz)

FetchContent_MakeAvailable(gtest)
