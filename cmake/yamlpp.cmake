include(FetchContent)
set(FETCHCONTENT_QUIET TRUE)

if(CANTERA_FOUND)
  return()
endif()

set(PACKAGE_NAME yaml-cpp)
set(REPO_URL "https://github.com/jbeder/yaml-cpp")
set(REPO_TAG "0.8.0")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} "" ON)

include_directories(${yaml-cpp_SOURCE_DIR}/include)
