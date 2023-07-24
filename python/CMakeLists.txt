# =============================================
# Wrap the C++ library and bind to Python
# =============================================

set(namel pyathena)
file(GLOB src_files
  pyathena.cpp
  )
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# global library
add_library(search_path OBJECT ${CMAKE_BINARY_DIR}/tests/globals.cpp)

# Set module targets

pybind11_add_module(${namel} MODULE ${src_files})

target_include_directories(
  ${namel}
  PRIVATE ${PYDISORT_INCLUDE_DIR} ${CMAKE_BINARY_DIR}
          ${PYTHON_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR})

target_link_libraries(${namel} PRIVATE
  search_path
  ${ATHENAPP_LIBRARY_${buildu}}
  ${CANOE_LIBRARY_${buildu}}
  ${NETCDF_LIBRARIES}
  ${PYDISORT_LIBRARY_${buildu}}
  ${PYTHON_LIBRARY_RELEASE}
  )

set_target_properties(${namel} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")

# Install configuration

install(
  TARGETS ${namel}
  COMPONENT python
  LIBRARY DESTINATION ${namel})
