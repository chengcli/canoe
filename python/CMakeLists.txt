# =========================================
# Wrap all C++ libraries and bind to Python
# =========================================

string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# global library
add_library(search_path OBJECT ${CMAKE_BINARY_DIR}/tests/globals.cpp)

# python site
get_filename_component(PYTHON_PARENT_DIR ${Python_SITELIB} DIRECTORY)
get_filename_component(PYTHON_VERSION_DIR ${PYTHON_PARENT_DIR} NAME)

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
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/site-packages)
endmacro()

# add python modules
add_python_module(pyathena)
add_python_module(pycanoe)
add_python_module(pyharp)
add_python_module(pyopacity)
add_python_module(utilities)
