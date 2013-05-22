cmake_minimum_required (VERSION 2.6)

FUNCTION(_Feel_COMPILER_DUMPVERSION _OUTPUT_VERSION)

  EXEC_PROGRAM(gcc
    ARGS  -dumpversion
    OUTPUT_VARIABLE _feel_COMPILER_VERSION
  )
  STRING(REGEX REPLACE "([0-9])\\.([0-9])(\\.[0-9])?" "\\1\\2"
    _feel_COMPILER_VERSION ${_feel_COMPILER_VERSION})

  SET(${_OUTPUT_VERSION} ${_feel_COMPILER_VERSION} PARENT_SCOPE)
ENDFUNCTION()


_Feel_COMPILER_DUMPVERSION(_feel_COMPILER_VERSION)

#MESSAGE(STATUS "feel_COMPILER_VERSION: gcc${_feel_COMPILER_VERSION}")
# MESSAGE(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")
#MESSAGE(STATUS "CMAKE_SYSTEM NAME: ${CMAKE_SYSTEM_NAME}")
#MESSAGE(STATUS "CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
#MESSAGE(STATUS "CMAKE_SYSTEM_VERSION: ${CMAKE_SYSTEM_VERSION}")
#set( build_name "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-gcc${_feel_COMPILER_VERSION}")
#MESSAGE(STATUS "build_name: ${build_name}")


find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)
set(CTEST_BUILD_NAME        "${osname}-${cpu}-gcc${_feel_COMPILER_VERSION}")


SET(MODEL Nightly)
IF(${CTEST_SCRIPT_ARG} MATCHES Experimental)
  SET(MODEL Experimental)
ENDIF()
IF(${CTEST_SCRIPT_ARG} MATCHES Continuous)
  SET(MODEL Continuous)
ENDIF()
MESSAGE( STATUS "Model: ${MODEL}" )

if ( ${MODEL} MATCHES Continuous )
  SET(FEELPP_ENABLE_ALL_DEFAULT OFF)
else()
  SET(FEELPP_ENABLE_ALL_DEFAULT ON)
endif()

SET (CTEST_INITIAL_CACHE "
// Enable tests
FEELPP_ENABLE_ALL:BOOL=${FEELPP_ENABLE_ALL_DEFAULT}
CMAKE_CXX_FLAGS:STRING=-std=c++0x -O3 -DOPTIMIZE -DNDEBUG -DNDEBUG_OLD
CMAKE_C_FLAGS:STRING=-std=c++0x -O3 -DOPTIMIZE -DNDEBUG -DNDEBUG_OLD
")


# -----------------------------------------------------------
# -- build specific
# -----------------------------------------------------------

## -- make command
## -----------------
find_program(MAKE NAMES make)

## -- Build options
# set(OPTION_BUILD                        "-j2")

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/sources/feel")
set(CTEST_BINARY_DIRECTORY  "$ENV{HOME}/sources/feel-${CTEST_BUILD_NAME}-${MODEL}")
set(CTEST_COMMAND "ctest -D ${MODEL}" )
SET(CTEST_CMAKE_COMMAND "cmake" )

SET (CTEST_GIT_COMMAND    "git" )
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  #set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone https://code.google.com/p/feelpp/ ${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone https://github.com/feelpp/feelpp.git ${CTEST_SOURCE_DIRECTORY}")
endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

#SET (CTEST_SVN_CHECKOUT   "${CTEST_SVN_COMMAND} co svn://scm.forge.imag.fr/var/lib/gforge/chroot/scmrepos/svn/life/trunk/life/trunk ${CTEST_SOURCE_DIRECTORY}")
#set (CTEST_UPDATE_COMMAND "${CTEST_SVN_COMMAND}")

# set(CTEST_BUILD_COMMAND     "make -j2")

# -----------------------------------------------------------
# -- commands
# -----------------------------------------------------------

set(CTEST_BUILD_COMMAND                "${MAKE} ${OPTION_BUILD}")

# -----------------------------------------------------------
# -- Settings
# -----------------------------------------------------------

set(CTEST_TIMEOUT           "600")

# default
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE 1)

if (${MODEL} MATCHES Nightly )

  # should ctest wipe the binary tree before running
  SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

endif()

if ( ${MODEL} MATCHES Continuous )
  # 10h duration
  while (${CTEST_ELAPSED_TIME} LESS 36000)
    # do work every 10 minutes if previous finished
    set (START_TIME ${CTEST_ELAPSED_TIME})
    ctest_start("Continuous")

    SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE 1)
#    ctest_update()
#    ctest_configure()
    # make sure that the target feel builds fully
#    ctest_build(TARGET feel)
#    ctest_submit()
    ctest_sleep( ${START_TIME} 300 ${CTEST_ELAPSED_TIME})
  endwhile()
endif()
