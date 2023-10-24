include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

option(PYTHON_BINDINGS "Build Python Bindings" OFF)

if(PYTHON_BINDINGS)
  # Execute python3-config to get the include flags
  find_package(Python COMPONENTS Interpreter Development)
  set(PYTHON_INCLUDE_DIR
      ${_Python_INCLUDE_DIR}
      CACHE PATH "include directory of python")
  set(PYTHON_EXECUTABLE
      ${_Python_EXECUTABLE}
      CACHE FILEPATH "executable of python")
  set(PYTHON_LIBRARY_RELEASE
      ${_Python_LIBRARY_RELEASE}
      CACHE FILEPATH "library of python")
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "import sysconfig; print(sysconfig.get_paths()['purelib'])"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  FetchContent_Declare(
    pybind11
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL https://github.com/pybind/pybind11/archive/v2.10.0.tar.gz)

  FetchContent_MakeAvailable(pybind11)

  include_directories(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()
