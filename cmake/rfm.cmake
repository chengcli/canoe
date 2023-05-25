include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)
set(GITHUB_USER $ENV{GH_ACCOUNT})
set(GITHUB_PAT $ENV{GH_TOKEN})

option(RFM "Build RFM" ON)

if(RFM)
  FetchContent_Declare(
    rfm
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL https://${GITHUB_USER}:${GITHUB_PAT}@github.com/luminoctum/rfm/archive/refs/tags/v4.33.tar.gz
  )

  FetchContent_GetProperties(rfm)

  if(NOT rfm_POPULATED)
    FetchContent_Populate(rfm)
    add_subdirectory(${rfm_SOURCE_DIR} ${rfm_BINARY_DIR})
  endif()
endif()
