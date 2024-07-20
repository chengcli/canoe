include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set(PACKAGE_NAME application)
set(REPO_URL "https://github.com/chengcli/application")
set(REPO_TAG "v0.7")
set(REPO_PATCH "None")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} ${REPO_PATCH} ON)
