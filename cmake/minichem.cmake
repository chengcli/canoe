include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set(PACKAGE_NAME minichem)
set(REPO_URL "https://github.com/chengcli/mini_chem")
set(REPO_TAG "3372a500f038f83228c9f8f944b3fb6b2dedc572")
set(REPO_PATCH "None")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} ${REPO_PATCH} OFF)
