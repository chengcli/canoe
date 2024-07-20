include(FetchContent)
set(FETCHCONTENT_QUIET TRUE)

if(CANTERA_FOUND)
  return()
endif()

set(PACKAGE_NAME yaml-cpp)
set(REPO_URL "https://github.com/jbeder/yaml-cpp")
set(REPO_TAG "0.8.0")
set(REPO_PATCH "None")
add_package(${PACKAGE_NAME} ${REPO_URL} ${REPO_TAG} ${REPO_PATCH} ON)

include_directories(${yaml-cpp_SOURCE_DIR}/include)
set(YAML_CPP_LIBRARIES yaml-cpp)
