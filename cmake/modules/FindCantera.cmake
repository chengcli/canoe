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

  set(SOURCE_PATH "${CANTERA_INCLUDE_DIR}/cantera/ext/")
  set(LINK_PATH "ext1")
  if(UNIX)
    execute_process(
      COMMAND ln -sf ${SOURCE_PATH} ${LINK_PATH}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      RESULT_VARIABLE result)
  elseif(WIN32)
    execute_process(
      COMMAND cmd.exe /c mklink /D ${LINK_PATH} ${SOURCE_PATH}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      RESULT_VARIABLE result)
  else()
    message(
      FATAL_ERROR "Symbolic link creation is not supported on this platform.")
  endif()

  if(result)
    message(FATAL_ERROR "Failed to create symbolic link: ${LINK_PATH}")
  else()
    message(STATUS "Symbolic link created: ${LINK_PATH} -> ${SOURCE_PATH}")
  endif()
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/ext1)
else()
  set(CANTERA_INCLUDE_DIR "")
  set(CANTERA_LIBRARY "")
  set(CANTERA_FOUND FALSE)
endif()
