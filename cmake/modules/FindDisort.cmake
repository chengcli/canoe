# Find the disort includes and library
#
# DISORT_INCLUDE_DIR  - where to find disort.hpp
#
# DISORT_LIBRARY      - Link these libraries when using DISORT
#
# DISORT_FOUND        - True if DISORT found
#
# Normal usage would be:
#
# find_package(Disort REQUIRED) include_directories(${DISORT_INCLUDE_DIR})
# target_link_libraries(${DISORT_LIBRARY})

include(FindPackageHandleStandardArgs)

macro(__disort_determine_version)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c
            "import pydisort; print(pydisort.__version__)"
    OUTPUT_VARIABLE DISORT_PEP440_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE result)

  if(result GREATER 0)
    message(FATAL_ERROR "cannot determine PEP 440 version of Pydisort!")
  endif()

  if(DISORT_PEP440_VERSION MATCHES "^[0-9]+\.[0-9]+(\.[0-9]+)?")
    set(DISORT_VERSION ${CMAKE_MATCH_0})
  endif()

  unset(result)
endmacro()

# Find Python
set(Python3_FIND_VIRTUALENV ONLY)
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import pydisort; print(pydisort.__file__)"
    OUTPUT_VARIABLE disort_init_file
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

  cmake_path(CONVERT ${disort_init_file} TO_CMAKE_PATH_LIST disort_init_file
             NORMALIZE)

  cmake_path(REPLACE_FILENAME disort_init_file lib OUTPUT_VARIABLE
             disort_lib_dir)
  cmake_path(REPLACE_FILENAME disort_init_file disort OUTPUT_VARIABLE
             disort_include_dir)

  unset(disort_init_file)
endif()

# Step 1: Find header
find_path(
  _DISORT_HEADER_DIR disort.hpp
  HINTS ${disort_include_dir}
        /opt/homebrew/include
        /usr/include
        DISORT_DIR/include
        DISORT_INC
        $ENV{DISORT_INC}
        $ENV{DISORT_DIR}/include
        $ENV{DISORT_ROOT}/include)

# Step 2: Go one level up from the header location
get_filename_component(DISORT_INCLUDE_DIR ${_DISORT_HEADER_DIR} DIRECTORY)

# Step 3: Save to cache
set(DISORT_INCLUDE_DIR
    "${DISORT_INCLUDE_DIR}"
    CACHE FILEPATH "Path to a file.")

# Step 4: unset the internal temp variable
unset(_DISORT_HEADER_DIR CACHE)

find_library(
  DISORT_LIBRARY disort_release
  HINTS ${disort_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        DISORT_DIR/lib
        DISORT_LIB
        $ENV{DISORT_LIB}
        $ENV{DISORT_DIR}/lib
        $ENV{DISORT_ROOT}/lib)

set(disort_required_vars DISORT_LIBRARY DISORT_INCLUDE_DIR)
mark_as_advanced(${disort_required_vars})

if(DISORT_LIBRARY)
  __disort_determine_version()
endif()

find_package_handle_standard_args(
  Disort
  REQUIRED_VARS ${disort_required_vars}
  VERSION_VAR DISORT_VERSION)

unset(disort_lib_dir)
unset(disort_include_dir)
unset(disort_required_vars)
