# Find the harp includes and library
#
# HARP_INCLUDE_DIR  - where to find element.hpp
#
# HARP_LIBRARY      - Link these libraries when using HARP
#
# HARP_FOUND        - True if harp found
#
# Normal usage would be:
#
# find_package(Harp REQUIRED) include_directories(${HARP_INCLUDE_DIR})
# target_link_libraries(${HARP_LIBRARY})

include(FindPackageHandleStandardArgs)

macro(__harp_determine_version)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import pyharp; print(pyharp.__version__)"
    OUTPUT_VARIABLE HARP_PEP440_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE result)

  if(result GREATER 0)
    message(FATAL_ERROR "cannot determine PEP 440 version of Pyharp!")
  endif()

  if(HARP_PEP440_VERSION MATCHES "^[0-9]+\.[0-9]+(\.[0-9]+)?")
    set(HARP_VERSION ${CMAKE_MATCH_0})
  endif()

  unset(result)
endmacro()

# Find Python
set(Python3_FIND_VIRTUALENV ONLY)
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import pyharp; print(pyharp.__file__)"
    OUTPUT_VARIABLE harp_init_file
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

  cmake_path(CONVERT ${harp_init_file} TO_CMAKE_PATH_LIST harp_init_file
             NORMALIZE)

  cmake_path(REPLACE_FILENAME harp_init_file lib OUTPUT_VARIABLE harp_lib_dir)
  cmake_path(REPLACE_FILENAME harp_init_file harp OUTPUT_VARIABLE
             harp_include_dir)

  unset(harp_init_file)
endif()

# Step 1: Find header
find_path(
  _HARP_HEADER_DIR element.hpp
  HINTS ${harp_include_dir}
        /opt/homebrew/include
        /usr/include
        HARP_DIR/include
        HARP_INC
        $ENV{HARP_INC}
        $ENV{HARP_DIR}/include
        $ENV{HARP_ROOT}/include)

# Step 2: Go one level up from the header location
get_filename_component(HARP_INCLUDE_DIR ${_HARP_HEADER_DIR} DIRECTORY)

# Step 3: Save to cache
set(HARP_INCLUDE_DIR
    "${HARP_INCLUDE_DIR}"
    CACHE FILEPATH "Path to a file.")

# Step 4: unset the internal temp variable
unset(_HARP_HEADER_DIR CACHE)

find_library(
  HARP_LIBRARY harp_release
  HINTS ${harp_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        HARP_DIR/lib
        HARP_LIB
        $ENV{HARP_LIB}
        $ENV{HARP_DIR}/lib
        $ENV{HARP_ROOT}/lib)

set(harp_required_vars HARP_LIBRARY HARP_INCLUDE_DIR)
mark_as_advanced(${harp_required_vars})

if(APPLE)
  link_directories(${Python3_SITELIB}/pyharp/.dylibs)
else()
  link_directories(${Python3_SITELIB}/pyharp.libs)
endif()

if(HARP_LIBRARY)
  __harp_determine_version()
endif()

find_package_handle_standard_args(
  Harp
  REQUIRED_VARS ${harp_required_vars}
  VERSION_VAR HARP_VERSION)

unset(harp_lib_dir)
unset(harp_include_dir)
unset(harp_required_vars)
