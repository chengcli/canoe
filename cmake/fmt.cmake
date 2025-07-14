include(FetchContent)
set(FETCHCONTENT_QUIET TRUE)

set(PACKAGE_NAME fmt)
set(REPO_URL "https://github.com/fmtlib/fmt")
set(REPO_TAG "11.1.2")

add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} "" ON)
set(FMT_INCLUDE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/_deps/${PACKAGE_NAME}-src/include"
    CACHE PATH "fmt include directory")
