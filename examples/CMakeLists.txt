# ========================
# A collection of examples
# ========================

if(${HYDROSTATIC})
  #add_subdirectory(1d-rad-jupiter)
  add_subdirectory(2023-Li-uranus)
  #add_subdirectory(2018-Li-harp)
else()
  if(${EQUATION_OF_STATE} STREQUAL "shallow_xy")
    add_subdirectory(2020-Li-polar-vortex)
  else()
    add_subdirectory(2024-Li-plume)
  endif()
endif()

if (${TASK} STREQUAL "juno")
  add_subdirectory(2023-jupiter-mwr-eq)
endif()

if (${TASK} STREQUAL "saturn_radio")
  add_subdirectory(2023-Li-saturn-vla)
endif()

if (${TASK} STREQUAL "straka" OR
    ${TASK} STREQUAL "robert" OR
    ${TASK} STREQUAL "bryan")
  add_subdirectory(2019-Li-snap)
endif()

if (${TASK} STREQUAL "ktable_jup")
  add_subdirectory(2024-FDing-jupiter-rt)
endif()

if (${TASK} STREQUAL "ktable_earth")
  add_subdirectory(2024-CLi-earth-rt)
endif()

if (${TASK} STREQUAL "exo2" OR
    ${TASK} STREQUAL "exo3")
  add_subdirectory(2023-Chen-exo3)
endif()

if (${TASK} STREQUAL "hjupiter")
  add_subdirectory(2024-XZhang-cloud-rt)
endif()

if (${TASK} STREQUAL "rcemip")
  add_subdirectory(2024-CLi-RCEMIP)
endif()
