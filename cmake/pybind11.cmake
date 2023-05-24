include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  pybind11
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/pybind/pybind11/archive/v2.10.0.tar.gz)

FetchContent_GetProperties(pybind11)

if(NOT pybind11_POPULATED)
  FetchContent_Populate(pybind11)
  add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()
