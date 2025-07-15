# Find the kintera includes and library
#
# KINTERA_INCLUDE_DIR:  Where to find reaction.hpp
#
# KINTERA_LIBRARY:      Link these libraries when using KINTERA
#
# VAPORS_LIBRARY:       Link these libraries when using KINTERA
#
# KINTERA_FOUND:        True if kintera found
#
# Normal usage would be:
#
# find_package(Kintera REQUIRED)
#
# include_directories(${KINTERA_INCLUDE_DIR})
#
# target_link_libraries(${KINTERA_LIBRARY} ${VAPORS_LIBRARY})

include(FindPackageHandleStandardArgs)

macro(__kintera_determine_version)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c
            "import kintera; print(kintera.__version__)"
    OUTPUT_VARIABLE KINTERA_PEP440_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE result)

  if(result GREATER 0)
    message(FATAL_ERROR "cannot determine PEP 440 version of Kintera!")
  endif()

  if(KINTERA_PEP440_VERSION MATCHES "^[0-9]+\.[0-9]+(\.[0-9]+)?")
    set(KINTERA_VERSION ${CMAKE_MATCH_0})
  endif()

  unset(result)
endmacro()

# Find Python
set(Python3_FIND_VIRTUALENV ONLY)
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import kintera; print(kintera.__file__)"
    OUTPUT_VARIABLE kintera_init_file
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

  cmake_path(CONVERT ${kintera_init_file} TO_CMAKE_PATH_LIST kintera_init_file
             NORMALIZE)

  cmake_path(REPLACE_FILENAME kintera_init_file lib OUTPUT_VARIABLE
             kintera_lib_dir)
  cmake_path(REPLACE_FILENAME kintera_init_file kintera OUTPUT_VARIABLE
             kintera_include_dir)

  unset(kintera_init_file)
endif()

# Step 1: Find header
find_path(
  _KINTERA_HEADER_DIR reaction.hpp
  HINTS ${kintera_include_dir}
        /opt/homebrew/include
        /usr/include
        KINTERA_DIR/include
        KINTERA_INC
        $ENV{KINTERA_INC}
        $ENV{KINTERA_DIR}/include
        $ENV{KINTERA_ROOT}/include)

# Step 2: Go one level up from the header location
get_filename_component(KINTERA_INCLUDE_DIR ${_KINTERA_HEADER_DIR} DIRECTORY)

# Step 3: Save to cache
set(KINTERA_INCLUDE_DIR
    "${KINTERA_INCLUDE_DIR}"
    CACHE FILEPATH "Path to a file.")

# Step 4: unset the internal temp variable
unset(_KINTERA_HEADER_DIR CACHE)

find_library(
  KINTERA_LIBRARY kintera_release
  HINTS ${kintera_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        KINTERA_DIR/lib
        KINTERA_LIB
        $ENV{KINTERA_LIB}
        $ENV{KINTERA_DIR}/lib
        $ENV{KINTERA_ROOT}/lib)

find_library(
  KINTERA_CUDA_LIBRARY kintera_cuda_release
  HINTS ${kintera_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        KINTERA_DIR/lib
        KINTERA_LIB
        $ENV{KINTERA_LIB}
        $ENV{KINTERA_DIR}/lib
        $ENV{KINTERA_ROOT}/lib)

find_library(
  VAPORS_LIBRARY vapors_release
  HINTS ${kintera_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        KINTERA_DIR/lib
        KINTERA_LIB
        $ENV{KINTERA_LIB}
        $ENV{KINTERA_DIR}/lib
        $ENV{KINTERA_ROOT}/lib)

if(${CUDAToolKit_FOUND})
  set(kintera_required_vars VAPORS_LIBRARY KINTERA_LIBRARY KINTERA_CUDA_LIBRARY
                            KINTERA_INCLUDE_DIR)
else()
  set(kintera_required_vars VAPORS_LIBRARY KINTERA_LIBRARY KINTERA_INCLUDE_DIR)
endif()

mark_as_advanced(${kintera_required_vars})

if(KINTERA_LIBRARY)
  __kintera_determine_version()
endif()

find_package_handle_standard_args(
  Kintera
  REQUIRED_VARS ${kintera_required_vars}
  VERSION_VAR KINTERA_VERSION)

unset(kintera_lib_dir)
unset(kintera_include_dir)
unset(kintera_required_vars)
