macro(add_package name url tag option)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(Yellow "${Esc}[33m")

  set(CACHE_DIR "${CMAKE_SOURCE_DIR}/.cache")
  set(CACHE_FILE "${CACHE_DIR}/${name}-${tag}.tar.gz")

  string(TOUPPER ${name} nameu)
  option(${nameu} "Build ${name}" ${option})

  if(${nameu})
    if(EXISTS ${CACHE_FILE})
      if(NOT EXISTS ${CMAKE_BINARY_DIR}/_deps/${name}-src)
        message(
          STATUS "${Yellow}Find cached file ${name}-${tag}.tar.gz${ColorReset}")
        execute_process(COMMAND tar -xzf ${CACHE_FILE} -C ${CACHE_DIR})
        execute_process(COMMAND mv ${CACHE_DIR}/${name}-src
                                ${CMAKE_BINARY_DIR}/_deps/)
      endif()

      FetchContent_Declare(${name} SOURCE_DIR
                                   ${CMAKE_BINARY_DIR}/_deps/${name}-src)
    else()
      FetchContent_Declare(
        ${name}
        GIT_REPOSITORY ${url}
        GIT_TAG ${tag})
    endif()

    FetchContent_MakeAvailable(${name})
    include_directories(${${name}_SOURCE_DIR})

    if(NOT EXISTS ${CACHE_FILE})
      message(STATUS "Creating ${CACHE_FILE}")
      execute_process(COMMAND tar -czf ${CACHE_DIR}/${name}-${tag}.tar.gz -C
                              ${CMAKE_BINARY_DIR}/_deps ${name}-src)
    endif()
  else()
    message(STATUS "${Yellow}Not building ${name}${ColorReset}")
  endif()
endmacro()
