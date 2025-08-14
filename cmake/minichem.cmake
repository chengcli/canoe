include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set(patch_command
    git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/35.mini_ch_i_dlsode.f90.patch)

set(PACKAGE_NAME minichem)
set(REPO_URL "https://github.com/chengcli/mini_chem")
set(REPO_TAG "3372a500f038f83228c9f8f944b3fb6b2dedc572")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} "${patch_command}" OFF)
