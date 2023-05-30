# ============================
# A collection of unit tests
# ============================

# Enable testing functionality
enable_testing()

if (${HYDROSTATIC})
  add_subdirectory(1d-rad-jupiter)
endif()
