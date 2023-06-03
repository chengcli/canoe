include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

FetchContent_Declare(
  gtest
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.tar.gz)

FetchContent_GetProperties(gtest)

if(NOT gtest_POPULATED)
  FetchContent_Populate(gtest)
  add_subdirectory(${gtest_SOURCE_DIR} ${gtest_BINARY_DIR})
endif()

FetchContent_MakeAvailable(gtest)
