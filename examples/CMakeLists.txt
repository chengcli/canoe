# ========================
# A collection of examples
# ========================

if (${TASK} STREQUAL "juno")
  add_subdirectory(2024-JHu-juno-mwr)
endif()

if (${TASK} STREQUAL "uranus_mwr")
  add_subdirectory(2023-Li-uranus)
endif()

if (${TASK} STREQUAL "saturn_radio")
  add_subdirectory(2023-Li-saturn-vla)
endif()

if (${TASK} STREQUAL "straka" OR
    ${TASK} STREQUAL "robert" OR
    ${TASK} STREQUAL "bryan" OR
    ${TASK} STREQUAL "jcrm")
  add_subdirectory(2019-Li-snap)
endif()

if (${TASK} STREQUAL "jupiter_dry" OR
    ${TASK} STREQUAL "jupiter_h2o" OR
    ${TASK} STREQUAL "jupiter_crm")
  add_subdirectory(2024-YZhou-jupiter)
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

if (${TASK} STREQUAL "rtdry" OR
    ${TASK} STREQUAL "rtcld" OR
    ${TASK} STREQUAL "rtcldexo")
  add_subdirectory(2024-XZhang-cloud-rt)
endif()

if (${TASK} STREQUAL "minichem")
  add_subdirectory(2024-XZhang-minichem)
endif()

if (${TASK} STREQUAL "rcemip")
  add_subdirectory(2024-CLi-RCEMIP)
endif()

if (${TASK} STREQUAL "vortex")
  add_subdirectory(2020-Li-polar-vortex)
endif()

if (${TASK} STREQUAL "amars")
  add_subdirectory(2024-CMetz-amars)
endif()

if (${TASK} STREQUAL "geyser")
  add_subdirectory(2024-YixiaoZ-geyser)
endif()
