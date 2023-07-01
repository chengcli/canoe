# * Try to find CFITSIO. Once executed, this module will define: Variables
#   defined by this module: CFITSIO_FOUND        - system has CFITSIO
#   CFITSIO_INCLUDE_DIR  - the CFITSIO include directory (cached)
#   CFITSIO_INCLUDE_DIRS - the CFITSIO include directories (identical to
#   CFITSIO_INCLUDE_DIR) CFITSIO_LIBRARY      - the CFITSIO library (cached)
#   CFITSIO_LIBRARIES    - the CFITSIO libraries (identical to CFITSIO_LIBRARY)
#
# This module will use the following enviornmental variable when searching for
# CFITSIO: CFITSIO_ROOT_DIR     - CFITSIO root directory
#

if(NOT cfitsio_FOUND)
  find_path(
    CFITSIO_INCLUDE_DIR fitsio.h
    HINTS $ENV{CFITSIO_ROOT_DIR} /usr /usr/local/
    PATH_SUFFIXES include include/cfitsio)
  find_library(
    CFITSIO_LIBRARY cfitsio /usr/lib /usr/lib64
    HINTS $ENV{CFITSIO_ROOT_DIR}
    PATH_SUFFIXES lib)

  mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(cfitsio DEFAULT_MSG CFITSIO_LIBRARY
                                    CFITSIO_INCLUDE_DIR)

  set(CFITSIO_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
  set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY})

endif(NOT cfitsio_FOUND)
