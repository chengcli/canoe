include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

FetchContent_Declare(
  yaml-cpp
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz)

FetchContent_GetProperties(yaml-cpp)

if(NOT yaml-cpp_POPULATED)
  message(STATUS "Fetching yaml-cpp...")
  FetchContent_Populate(yaml-cpp)
  add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
endif()

# FetchContent_MakeAvailable(yamlcpp)

# Where yaml-cpp's .h files can be found.
include_directories(${yaml-cpp_SOURCE_DIR}/include)
