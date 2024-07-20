# Fetch disort and build
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set_if_empty(ACCOUNT $ENV{GH_ACCOUNT})
set_if_empty(TOKEN $ENV{GH_TOKEN})

set(PACKAGE_NAME disort)
set(REPO_URL "https://${ACCOUNT}:${TOKEN}@github.com/zoeyzyhu/pydisort")
set(REPO_TAG "v0.8.0")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} "" OFF)
