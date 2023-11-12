# =========================================
# Wrap all C++ libraries and bind to Python
# =========================================

string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# global library
add_library(search_path OBJECT ${CMAKE_BINARY_DIR}/tests/globals.cpp)

# Execute a Python script to get the site-packages path
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

macro(add_python_module namel)
  # set module targets
  pybind11_add_module(${namel} MODULE ${namel}.cpp)

  target_include_directories(${namel}
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CPPDISORT_INCLUDE_DIR}
    ${CANOE_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIR}
    ${PYBIND11_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    )

  target_link_libraries(${namel} PRIVATE
    search_path
    ${CANOE_LIBRARY_${buildu}}
    )

  set_target_properties(${namel}
    PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}")

  # Install configuration
  install(
    TARGETS ${namel}
    COMPONENT python
    DESTINATION ${PYTHON_SITE_PACKAGES})
endmacro()

# add python modules
add_python_module(pyathena)
add_python_module(pycanoe)
add_python_module(pyharp)
add_python_module(utilities)
