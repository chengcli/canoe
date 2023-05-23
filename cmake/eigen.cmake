include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

FetchContent_Declare(
  eigen URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE)

FetchContent_GetProperties(eigen)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)

if(NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
  add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR})
endif()

set(EIGEN_INCLUDE_DIR
    ${Eigen3_SOURCE_DIR}
    CACHE PATH "Eigen include directory")

# find_package(Eigen3 3.3 REQUIRED NO_MODULE) include(FetchContent)
# FetchContent_Declare( Eigen GIT_REPOSITORY
# https://gitlab.com/libeigen/eigen.git GIT_TAG master GIT_SHALLOW TRUE
# GIT_PROGRESS TRUE)

# set(CMAKE_POLICY_DEFAULT_CMP0077 NEW) set(EIGEN_BUILD_DOC OFF)
# set(EIGEN_BUILD_PKGCONFIG OFF) FetchContent_MakeAvailable(Eigen)
