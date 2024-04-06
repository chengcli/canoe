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
    ${CMAKE_CURRENT_SOURCE_DIR}/patches/27.coordinates.patch)

# Define the path to the local directory
set(local_dir "$ENV{HOME}/modules/application")

# Check if the local directory exists
if(EXISTS ${local_dir})
    message(STATUS "Using local directory: ${local_dir}")
    # Use FetchContent with the local directory
    FetchContent_Declare(
        application
        URL "file://${local_dir}"
        URL_HASH "SHA256=0"  # Dummy hash to satisfy FetchContent
        PATCH_COMMAND ${patch_command}
        DOWNLOAD_COMMAND rsync -a ${local_dir}/ <SOURCE_DIR>
        UPDATE_DISCONNECTED TRUE)
else()
    message(STATUS "Local directory not found, downloading from GitHub")
    FetchContent_Declare(
        application
        GIT_REPOSITORY https://github.com/chengcli/application/archive/refs/tags/v0.7.tar.gz
        GIT_TAG snap-mods
        PATCH_COMMAND ${patch_command}
        UPDATE_DISCONNECTED TRUE)
# DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL
# https://github.com/chengcli/application/archive/refs/tags/v0.7.tar.gz)
endif()

FetchContent_MakeAvailable(application)

include_directories(${application_SOURCE_DIR})
