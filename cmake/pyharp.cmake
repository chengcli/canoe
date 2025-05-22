include(FetchContent)
set(FETCHCONTENT_QUIET TRUE)

set(PACKAGE_NAME pyharp)
set(REPO_URL "https://github.com/chengcli/pyharp")
set(REPO_TAG "v1.6.3")

option(BUILD_TESTS "Build tests" OFF)
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} "" ON)
include_directories(${pyharp_SOURCE_DIR})
