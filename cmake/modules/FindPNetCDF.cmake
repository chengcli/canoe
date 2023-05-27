# * Find PNetCDF Find the PNetcdf include and library
#
# PNETCDF_INCLUDE_DIR  - where to find pnetcdf.h PNETCDF_LIBRARY      - link
# library PNETCDF_FOUND        - True if PNETCDF found include requried
# interfaces
#
# Normal usage would be:

find_path(PNETCDF_INCLUDE_DIR pnetcdf.h HINTS $ENV{HOME}/opt/include)

find_library(
  PNETCDF_LIBRARY
  NAMES pnetcdf
  HINTS $ENV{HOME}/opt/lib)

if(PNETCDF_INCLUDE_DIR)
  mark_as_advanced(PNETCDF_LIBRARY)
  mark_as_advanced(PNETCDF_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PNetCDF DEFAULT_MSG PNETCDF_LIBRARY
                                    PNETCDF_INCLUDE_DIR)
endif()
