# ========================
# A collection of examples
# ========================

if(${HYDROSTATIC})
  #add_subdirectory(1d-rad-jupiter)
  add_subdirectory(2023-jupiter-mwr-eq)
  add_subdirectory(2023-Li-saturn-vla)
  add_subdirectory(2023-Li-uranus)
  #add_subdirectory(2018-Li-harp)
  if(${EQUATION_OF_STATE} STREQUAL "shallow_yz")
    add_subdirectory(2023-Chen-exo3)
  endif()
else()
  if(${EQUATION_OF_STATE} STREQUAL "shallow_xy")
    add_subdirectory(2020-Li-polar-vortex)
  else()
    add_subdirectory(2019-Li-snap)
    add_subdirectory(2024-Li-plume)
    add_subdirectory(2024-XZhang-cloud-rt)
  endif()
endif()
