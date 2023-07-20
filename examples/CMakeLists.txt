# ========================
# A collection of examples
# ========================

if(${HYDROSTATIC})
  #add_subdirectory(1d-rad-jupiter)
  add_subdirectory(2023-jupiter-mwr-eq)
else()
  if(${EQUATION_OF_STATE} STREQUAL "shallow_xy")
    add_subdirectory(2020-Li-polar-vortex)
  else()
    add_subdirectory(2019-Li-snap)
  endif()
endif()
