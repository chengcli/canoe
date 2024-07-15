# =========================================
# Wrap all C++ libraries and bind to Python
# =========================================
if (NOT PYTHON_BINDINGS)
  return()
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# global library
add_library(search_path OBJECT ${CMAKE_BINARY_DIR}/tests/globals.cpp)

# Execute a Python script to get the site-packages path
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(namel canoe)
set_if_empty(PGEN empty)

# set module targets
pybind11_add_module(${namel} MODULE
  pyathena.cpp
  pycanoe.cpp
  pyharp.cpp
  pysnap.cpp
  pyutils.cpp
  pyinversion.cpp
  ${PGEN}.cpp
  )

target_include_directories(${namel}
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CPPDISORT_INCLUDE_DIR}
  ${CANOE_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIR}
  ${PYBIND11_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_PATH}
  ${MPI_CXX_HEADER_DIR}
  )

target_link_libraries(${namel} PRIVATE
  search_path
  ${CANOE_LIBRARY_${buildu}}
  )

set_target_properties(${namel}
  PROPERTIES
  PREFIX "${PYTHON_MODULE_PREFIX}"
  SUFFIX "${PYTHON_MODULE_EXTENSION}")

# softlink library
execute_process(
  COMMAND ln -sf ${CMAKE_BINARY_DIR}/python/${namel}${PYTHON_MODULE_EXTENSION}
  ${CMAKE_BINARY_DIR}/lib/${namel}${PYTHON_MODULE_EXTENSION}
  )
