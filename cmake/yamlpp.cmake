include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

FetchContent_Declare(
  yamlcpp
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.7.0.tar.gz
)

FetchContent_MakeAvailable(yamlcpp)

# Where yaml-cpp's .h files can be found.
include_directories(${yamlcpp_SOURCE_DIR}/include)
