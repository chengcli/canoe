# Fetch disort and build
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

set_if_empty(ACCOUNT $ENV{GH_ACCOUNT})
set_if_empty(TOKEN $ENV{GH_TOKEN})

option(DISORT "Build DISORT" OFF)

if(DISORT)
  FetchContent_Declare(
    pydisort
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL https://${ACCOUNT}:${TOKEN}@github.com/zoeyzyhu/pydisort/archive/refs/tags/v0.8.0.tar.gz
  )

  FetchContent_MakeAvailable(pydisort)
endif()
