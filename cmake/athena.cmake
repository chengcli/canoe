include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

set(patch_command
    git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/19.calculate_fluxes.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/21.new_blockdt.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/23.meshblock.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/24.time_integrator.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/25.constant_acc.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/26.scalars_flux.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/26.adiabatic_hydro.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/27.coordinates.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/28.bvals_cc.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/29.mesh.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/30.bvals_base.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/31.task_list.patch
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/32.scalars.patch)

FetchContent_Declare(
  athenapp
  GIT_REPOSITORY https://github.com/chengcli/athenapp/
  GIT_TAG snap-mods
  PATCH_COMMAND ${patch_command}
  UPDATE_DISCONNECTED TRUE)
# DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL
# https://github.com/chengcli/athenapp/archive/refs/tags/v0.8.tar.gz)

FetchContent_MakeAvailable(athenapp)

include_directories(${athenapp_SOURCE_DIR})
