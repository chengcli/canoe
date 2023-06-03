include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

FetchContent_Declare(
  yamlpp
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.7.0.tar.gz
)

FetchContent_GetProperties(yamlpp)

if(NOT yamlpp_POPULATED)
  FetchContent_Populate(yamlpp)
  add_subdirectory(${yamlpp_SOURCE_DIR} ${yamlpp_BINARY_DIR})
endif()

# set(YAMLPP_INCLUDE_DIR ${yamlpp_SOURCE_DIR} CACHE PATH "yaml-cpp include
# directory")
