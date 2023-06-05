include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

FetchContent_Declare(
  eigen
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz)

FetchContent_GetProperties(eigen)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)

FetchContent_MakeAvailable(eigen)

include_directories(${Eigen3_SOURCE_DIR})
