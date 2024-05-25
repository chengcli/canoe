include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

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
set(local_dir "$ENV{HOME}/modules/yamlpp")

# Check if the local directory exists
if(EXISTS ${local_dir})
    message(STATUS "Using local directory: ${local_dir}")
    # Use FetchContent with the local directory
    FetchContent_Declare(
        yaml-cpp
        URL "file://${local_dir}"
        URL_HASH "SHA256=0"  # Dummy hash to satisfy FetchContent
        PATCH_COMMAND ${patch_command}
        DOWNLOAD_COMMAND rsync -a ${local_dir}/ <SOURCE_DIR>
        UPDATE_DISCONNECTED TRUE)
else()
    message(STATUS "Local directory not found, downloading from GitHub")
    FetchContent_Declare(
        yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz
        GIT_TAG snap-mods
        PATCH_COMMAND ${patch_command}
        UPDATE_DISCONNECTED TRUE)
# DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL
# https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz)
endif()

FetchContent_GetProperties(yaml-cpp)

if(NOT yaml-cpp_POPULATED)
  message(STATUS "Fetching yaml-cpp...")
  FetchContent_Populate(yaml-cpp)
  add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
endif()

# FetchContent_MakeAvailable(yamlcpp)

# Where yaml-cpp's .h files can be found.
include_directories(${yaml-cpp_SOURCE_DIR}/include)
