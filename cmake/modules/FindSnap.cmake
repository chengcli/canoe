# Find the snap includes and library
#
# SNAP_INCLUDE_DIR:  Where to find snap.h
#
# SNAP_LIBRARY:      Link these libraries when using SNAP
#
# BC_LIBRARY:         Link these libraries when using SNAP
#
# SNAP_FOUND:        True if snap found
#
# Normal usage would be:
#
# find_package(Snap REQUIRED)
#
# include_directories(${SNAP_INCLUDE_DIR})
#
# target_link_libraries(${SNAP_LIBRARY} ${BC_LIBRARY})

include(FindPackageHandleStandardArgs)

macro(__snap_determine_version)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import snapy; print(snapy.__version__)"
    OUTPUT_VARIABLE SNAP_PEP440_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE result)

  if(result GREATER 0)
    message(FATAL_ERROR "cannot determine PEP 440 version of Snap!")
  endif()

  if(SNAP_PEP440_VERSION MATCHES "^[0-9]+\.[0-9]+(\.[0-9]+)?")
    set(SNAP_VERSION ${CMAKE_MATCH_0})
  endif()

  unset(result)
endmacro()

# Find Python
set(Python3_FIND_VIRTUALENV ONLY)
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import snapy; print(snapy.__file__)"
    OUTPUT_VARIABLE snap_init_file
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

  cmake_path(CONVERT ${snap_init_file} TO_CMAKE_PATH_LIST snap_init_file
             NORMALIZE)

  cmake_path(REPLACE_FILENAME snap_init_file lib OUTPUT_VARIABLE snap_lib_dir)
  cmake_path(REPLACE_FILENAME snap_init_file snap OUTPUT_VARIABLE
             snap_include_dir)

  unset(snap_init_file)
endif()

# Step 1: Find header
find_path(
  _SNAP_HEADER_DIR snap.h
  HINTS ${snap_include_dir}
        /opt/homebrew/include
        /usr/include
        SNAP_DIR/include
        SNAP_INC
        $ENV{SNAP_INC}
        $ENV{SNAP_DIR}/include
        $ENV{SNAP_ROOT}/include)

# Step 2: Go one level up from the header location
get_filename_component(SNAP_INCLUDE_DIR ${_SNAP_HEADER_DIR} DIRECTORY)

# Step 3: Save to cache
set(SNAP_INCLUDE_DIR
    "${SNAP_INCLUDE_DIR}"
    CACHE FILEPATH "Path to a file.")

# Step 4: unset the internal temp variable
unset(_SNAP_HEADER_DIR CACHE)

find_library(
  SNAP_LIBRARY snap_release
  HINTS ${snap_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        SNAP_DIR/lib
        SNAP_LIB
        $ENV{SNAP_LIB}
        $ENV{SNAP_DIR}/lib
        $ENV{SNAP_ROOT}/lib)

find_library(
  SNAP_CUDA_LIBRARY snap_cuda_release
  HINTS ${snap_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        SNAP_DIR/lib
        SNAP_LIB
        $ENV{SNAP_LIB}
        $ENV{SNAP_DIR}/lib
        $ENV{SNAP_ROOT}/lib)

find_library(
  BC_LIBRARY bc_release
  HINTS ${snap_lib_dir}
        /opt/homebrew/lib
        /usr/lib/x86_64-linux-gnu/
        SNAP_DIR/lib
        SNAP_LIB
        $ENV{SNAP_LIB}
        $ENV{SNAP_DIR}/lib
        $ENV{SNAP_ROOT}/lib)

if(${CUDAToolKit_FOUND})
  set(snap_required_vars BC_LIBRARY SNAP_LIBRARY SNAP_CUDA_LIBRARY
                         SNAP_INCLUDE_DIR)
else()
  set(snap_required_vars BC_LIBRARY SNAP_LIBRARY SNAP_INCLUDE_DIR)
endif()

mark_as_advanced(${snap_required_vars})

if(SNAP_LIBRARY)
  __snap_determine_version()
endif()

find_package_handle_standard_args(
  Snap
  REQUIRED_VARS ${snap_required_vars}
  VERSION_VAR SNAP_VERSION)

unset(snap_lib_dir)
unset(snap_include_dir)
unset(snap_required_vars)
