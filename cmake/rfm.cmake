include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set_if_empty(ACCOUNT $ENV{GH_ACCOUNT})
set_if_empty(TOKEN $ENV{GH_TOKEN})

set(PACKAGE_NAME rfm)
set(REPO_URL "https://${ACCOUNT}:${TOKEN}@github.com/luminoctum/rfm")
set(REPO_TAG "v5.20.2")
set(REPO_PATCH "None")
add_package_noinclude(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} ${REPO_PATCH} OFF)
