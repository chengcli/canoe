include(FetchContent)
set(FETCHCONTENT_QUIET TRUE)

set(PACKAGE_NAME gtest)
set(REPO_URL "https://github.com/chengcli/gtest")
set(REPO_TAG "v1.13.0")
set(REPO_PATCH "None")
set(INSTALL_GTEST OFF)

add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} ${REPO_PATCH} ON)
