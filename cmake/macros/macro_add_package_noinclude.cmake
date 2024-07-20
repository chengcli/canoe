macro(add_package_noinclude name url tag patch option)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(Yellow "${Esc}[33m")

  set(CACHE_DIR "${CMAKE_SOURCE_DIR}/.cache")
  set(CACHE_FILE "${CACHE_DIR}/${name}-${tag}.tar.gz")

  string(TOUPPER ${name} nameu)
  option(${nameu} "Build ${name}" ${option})

  if(${${nameu}})
    if(EXISTS ${CACHE_FILE})
      if(NOT EXISTS ${CMAKE_BINARY_DIR}/_deps)
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps)
      endif()

      if(NOT EXISTS ${CMAKE_BINARY_DIR}/_deps/${name}-src)
        message(
          STATUS
            "${Yellow}Found cached library ${name}-${tag}.tar.gz${ColorReset}")
        execute_process(COMMAND tar -xzf ${CACHE_FILE} -C
                                ${CMAKE_BINARY_DIR}/_deps)
      endif()

      FetchContent_Declare(${name} SOURCE_DIR
                                   ${CMAKE_BINARY_DIR}/_deps/${name}-src)
    else()
      if(${patch} STREQUAL "None")
        FetchContent_Declare(
          ${name}
          GIT_REPOSITORY ${url}
          GIT_TAG ${tag})
      else()
        FetchContent_Declare(
          ${name}
          GIT_REPOSITORY ${url}
          GIT_TAG ${tag}
          PATCH_COMMAND ${patch}
          UPDATE_DISCONNECTED TRUE)
      endif()
    endif()

    FetchContent_MakeAvailable(${name})

    if(NOT EXISTS ${CACHE_FILE})
      message(
        STATUS
          "${Yellow}Creating cached library ${name}-${tag}.tar.gz${ColorReset}")
      execute_process(COMMAND tar -czf ${CACHE_DIR}/${name}-${tag}.tar.gz -C
                              ${CMAKE_BINARY_DIR}/_deps ${name}-src)
    endif()
  else()
    message(STATUS "${Yellow}Not building ${name}${ColorReset}")
  endif()
endmacro()
