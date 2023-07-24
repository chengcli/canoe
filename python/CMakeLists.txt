# =========================================
# Wrap all C++ libraries and bind to Python
# =========================================

string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# global library
add_library(search_path OBJECT ${CMAKE_BINARY_DIR}/tests/globals.cpp)

macro(add_python_module namel)
  # set module targets
  pybind11_add_module(${namel} MODULE ${namel}.cpp)

  target_include_directories(
    ${namel}
    PRIVATE ${PYDISORT_INCLUDE_DIR} ${CMAKE_BINARY_DIR} ${CANOE_INCLUDE_DIR}
            ${PYTHON_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR})

  target_link_libraries(${namel} PRIVATE
    search_path
    ${ATHENAPP_LIBRARY_${buildu}}
    ${CANOE_LIBRARY_${buildu}}
    ${NETCDF_LIBRARIES}
    ${PYDISORT_LIBRARY_${buildu}}
    ${PYTHON_LIBRARY_RELEASE}
    )

  set_target_properties(${namel}
    PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}")
endmacro()

# add python modules
add_python_module(pyathena)
add_python_module(pyharp)

# Install configuration

install(
  TARGETS ${namel}
  COMPONENT python
  LIBRARY DESTINATION ${namel})
