include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
    athenapp
    URL https://github.com/chengcli/athenapp/archive/refs/tags/??v21.0-dev.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

FetchContent_GetProperties(athenapp)

if(NOT athenapp_POPULATED)
    FetchContent_Populate(athenapp)
    add_subdirectory(${athenapp_SOURCE_DIR})
endif()
