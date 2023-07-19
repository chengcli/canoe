include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

set_if_empty(ACCOUNT $ENV{GH_ACCOUNT})
set_if_empty(TOKEN $ENV{GH_TOKEN})

option(RFM "Build RFM" OFF)

if(RFM)
  FetchContent_Declare(
    rfm
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL https://${ACCOUNT}:${TOKEN}@github.com/luminoctum/rfm/archive/refs/tags/v5.20.2.tar.gz
  )

  FetchContent_MakeAvailable(rfm)
endif()
