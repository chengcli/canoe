include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set(PACKAGE_NAME exofmsrt)
set(REPO_URL "https://github.com/chengcli/Exo-FMS_column_ck")
set(REPO_TAG "36ba9fad0f87339c3c7b83dd9bf2b53eaba8223f")
set(REPO_PATCH "None")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} ${REPO_PATCH} OFF)
