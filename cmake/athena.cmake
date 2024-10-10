include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set(patch_command
    git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/19.calculate_fluxes.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/21.new_blockdt.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/23.meshblock.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/24.time_integrator.patch
    # ${CMAKE_CURRENT_SOURCE_DIR}/patches/25.constant_acc.patch
    # ${CMAKE_CURRENT_SOURCE_DIR}/patches/26.scalars_flux.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/26.adiabatic_hydro.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/27.coordinates.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/29.mesh.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/30.bvals_base.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/31.task_list.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/32.scalars.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/33.bvals_cc_cpp.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/34.bvals_cc_hpp.patch)

set(PACKAGE_NAME athenapp)
set(REPO_URL "https://github.com/chengcli/athenapp")
# set(REPO_TAG "12b11a6004827ba4732c248c021be55b2906ab82")
set(REPO_TAG "1750ec93a0969a9d572d8e0d6a094583b7007bb0")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} "${patch_command}" ON)
