# * Find Cantera Find the Cantera include and library
#
# Define the following variables
#
# CANTERA_INCLUDE_DIR  - where to find cantera/kinetics.h CANTERA_LIBRARY      -
# link library CANTERA_FOUND        - True if Cantera found
#
# Normal usage would be
#
# find_package (Cantera REQUIRED) target_include_directories
# (${CANTERA_INCLUDE_DIR}) target_link_libraries (${CANTERA_LIBRARY})

find_path(CANTERA_INCLUDE_DIR cantera/kinetics.h HINTS $ENV{HOME}/opt/include
                                                       $ENV{CANTERA_DIR})

find_library(
  CANTERA_LIBRARY
  NAMES cantera
  HINTS $ENV{HOME}/opt/lib)

if(CANTERA_INCLUDE_DIR)
  mark_as_advanced(CANTERA_LIBRARY)
  mark_as_advanced(CANTERA_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Cantera DEFAULT_MSG CANTERA_LIBRARY
                                    CANTERA_INCLUDE_DIR)
endif()
