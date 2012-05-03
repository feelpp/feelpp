
####################################################################
#
# Usage:
#  - create a new folder, let's call it cdash
#  - in that folder, do:
#    ctest -S path/to/feel/cmake/dashboard/testsuite.cmake[,option1=value1[,option2=value2]]
#
# Options:
#  - FEELPP_CXX: compiler, eg.: g++-4.2
#      default: default c++ compiler
#  - FEELPP_SITE: eg, sd-25154, or the name of the contributor, etc.
#      default: hostname
#  - FEELPP_BUILD_STRING: a string which identify the system/compiler. It should be formed like that:
#        <OS_name>-<OS_version>-<arch>-<compiler-version>
#      with:
#        <OS_name> = opensuse, debian, osx, windows, cygwin, freebsd, solaris, etc.
#        <OS_version> = 11.1, XP, vista, leopard, etc.
#        <arch> = i386, x86_64, ia64, powerpc, etc.
#        <compiler-version> = gcc-4.3.2, icc-11.0, MSVC-2008, etc.
#  - FEELPP_EXPLICIT_VECTORIZATION: novec, SSE2, Altivec
#       default: SSE2 for x86_64 systems, novec otherwise
#       Its value is automatically appended to FEELPP_BUILD_STRING
#  - FEELPP_CMAKE_DIR: path to cmake executable
#  - FEELPP_MODE: dashboard model, can be Experimental, Nightly, or Continuous
#      default: Nightly
#  - FEELPP_WORK_DIR: directory used to download the source files and make the builds
#      default: folder which contains this script
#  - FEELPP_CMAKE_ARGS: additional arguments passed to cmake
#  - FEELPP_GENERATOR_TYPE: allows to overwrite the generator type
#      default: nmake (windows
#      See http://www.cmake.org/cmake/help/cmake2.6docs.html#section_Generators for a complete
#      list of supported generators.
#  - FEELPP_NO_UPDATE: allows to submit dash boards from local repositories
#      This might be interesting in case you want to submit dashboards
#      including local changes.
#  - CTEST_SOURCE_DIRECTORY: path to feel's feel (use a new and empty folder, not the one you are working on)
#      default: <FEELPP_WORK_DIR>/feel
#  - CTEST_BINARY_DIRECTORY: build directory
#      default: <FEELPP_WORK_DIR>/nightly-<FEELPP_CXX>
#
# Here is an example running several compilers on a linux system:
# #!/bin/bash
# ARCH=`uname -m`
# SITE=`hostname`
# VERSION=opensuse-11.1
# WORK_DIR=/home/prudhomm/scratch/cdash
# # get the last version of the script
# wget http://bitbucket.org/feel/feel/raw/tip/test/testsuite.cmake -o $WORK_DIR/testsuite.cmake
# COMMON="ctest -S $WORK_DIR/testsuite.cmake,FEELPP_WORK_DIR=$WORK_DIR,FEELPP_SITE=$SITE,FEELPP_MODE=$1,FEELPP_BUILD_STRING=$OS_VERSION-$ARCH"
# $COMMON-gcc-3.4.6,FEELPP_CXX=g++-3.4
# $COMMON-gcc-4.0.1,FEELPP_CXX=g++-4.0.1
# $COMMON-gcc-4.3.2,FEELPP_CXX=g++-4.3,FEELPP_EXPLICIT_VECTORIZATION=novec
# $COMMON-gcc-4.3.2,FEELPP_CXX=g++-4.3,FEELPP_EXPLICIT_VECTORIZATION=SSE2
# $COMMON-icc-11.0,FEELPP_CXX=icpc
#
####################################################################

# process the arguments

set(ARGLIST ${CTEST_SCRIPT_ARG})
while(${ARGLIST} MATCHES  ".+.*")

  # pick first
  string(REGEX MATCH "([^,]*)(,.*)?" DUMMY ${ARGLIST})
  SET(TOP ${CMAKE_MATCH_1})

  # remove first
  string(REGEX MATCHALL "[^,]*,(.*)" DUMMY ${ARGLIST})
  SET(ARGLIST ${CMAKE_MATCH_1})

  # decompose as a pair key=value
  string(REGEX MATCH "([^=]*)(=.*)?" DUMMY ${TOP})
  SET(KEY ${CMAKE_MATCH_1})

  string(REGEX MATCH "[^=]*=(.*)" DUMMY ${TOP})
  SET(VALUE ${CMAKE_MATCH_1})

  # set the variable to the specified value
  if(VALUE)
    SET(${KEY} ${VALUE})
  else(VALUE)
    SET(${KEY} ON)
  endif(VALUE)

endwhile(${ARGLIST} MATCHES ".+.*")

####################################################################
# Automatically set some user variables if they have not been defined manually
####################################################################
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

if ( FEELPP_CTEST_CONFIG )
  include(${FEELPP_CTEST_CONFIG})
  message(STATUS "FEELPP_SITE: ${FEELPP_SITE}")
  message(STATUS "FEELPP_CXX: ${FEELPP_CXX}")
endif()

if(NOT FEELPP_SITE)
  site_name(FEELPP_SITE)
endif(NOT FEELPP_SITE)

if(NOT FEELPP_CMAKE_DIR)
  SET(FEELPP_CMAKE_DIR "")
endif(NOT FEELPP_CMAKE_DIR)

if (NOT FEELPP_CXX)
  set(FEELPP_CXX "g++")
endif(NOT FEELPP_CXX)

get_filename_component( FEELPP_CXX_NAME ${FEELPP_CXX} NAME )

if(NOT FEELPP_BUILD_STRING)

  # let's try to find all information we need to make the build string ourself

  # OS
  build_name(FEELPP_OS_VERSION)

  # arch
  set(FEELPP_ARCH ${CMAKE_SYSTEM_PROCESSOR})
  if(WIN32)
    set(FEELPP_ARCH $ENV{PROCESSOR_ARCHITECTURE})
  else(WIN32)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE FEELPP_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif(WIN32)


  set(FEELPP_BUILD_STRING ${FEELPP_OS_VERSION}${FEELPP_ARCH}-${FEELPP_CXX_NAME})

endif(NOT FEELPP_BUILD_STRING)

if(DEFINED FEELPP_EXPLICIT_VECTORIZATION)
  set(FEELPP_BUILD_STRING ${FEELPP_BUILD_STRING}-${FEELPP_EXPLICIT_VECTORIZATION})
endif(DEFINED FEELPP_EXPLICIT_VECTORIZATION)

if(NOT FEELPP_WORK_DIR)
  set(FEELPP_WORK_DIR ${CTEST_SCRIPT_DIRECTORY})
endif(NOT FEELPP_WORK_DIR)

if(NOT CTEST_SOURCE_DIRECTORY)
  SET (CTEST_SOURCE_DIRECTORY "${FEELPP_WORK_DIR}/feelpp")
endif(NOT CTEST_SOURCE_DIRECTORY)

if(NOT CTEST_BINARY_DIRECTORY)
  SET (CTEST_BINARY_DIRECTORY "${FEELPP_WORK_DIR}/${FEELPP_MODE}_${FEELPP_CXX_NAME}")
endif(NOT CTEST_BINARY_DIRECTORY)

if(NOT FEELPP_MODE)
  set(FEELPP_MODE Nightly)
endif(NOT FEELPP_MODE)

## mandatory variables (the default should be ok in most cases):

#if(NOT FEELPP_NO_UPDATE)
find_program(CTEST_GIT_COMMAND NAMES git)
find_program(CTEST_SVN_COMMAND NAMES svn)
#SET (CTEST_GIT_COMMAND "git")
#SET (CTEST_SVN_CHECKOUT   "${CTEST_GIT_COMMAND} co svn://scm.forge.imag.fr/var/lib/gforge/chroot/scmrepos/svn/life/trunk/life/trunk ${CTEST_SOURCE_DIRECTORY}")
#SET (CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone https://code.google.com/p/feelpp/")
set (CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

  #SET(CTEST_BACKUP_AND_RESTORE TRUE) # the backup is SVN related ...
#endif(NOT FEELPP_NO_UPDATE)
foreach(module ${FEELPP_MODULES})
  # update the modules using svn update
  execute_process(
    COMMAND "cd ${CTEST_SOURCE_DIRECTORY}/${module} && ${CTEST_SVN_COMMAND} update"
    OUTPUT_VARIABLE MODULE_OUTPUT)
  message(STATUS "updated ${module} : ${MODULE_OUTPUT}")
endforeach()
####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this make sure we get consistent outputs
SET($ENV{LC_MESSAGES} "en_EN")

if (UNIX)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif(UNIX)

# if(DEFINED FEELPP_EXPLICIT_VECTORIZATION)
#   if(FEELPP_EXPLICIT_VECTORIZATION MATCHES SSE2)
#     set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEELPP_TEST_SSE2=ON")
#   elseif(FEELPP_EXPLICIT_VECTORIZATION MATCHES SSE3)
#     set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEELPP_TEST_SSE2=ON -DFEELPP_TEST_SSE3=ON")
#   elseif(FEELPP_EXPLICIT_VECTORIZATION MATCHES SSSE3)
#     set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEELPP_TEST_SSE2=ON -DFEELPP_TEST_SSE3=ON -DFEELPP_TEST_SSSE3=ON")
#   elseif(FEELPP_EXPLICIT_VECTORIZATION MATCHES SSE4_1)
#     set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEELPP_TEST_SSE2=ON -DFEELPP_TEST_SSE3=ON -DFEELPP_TEST_SSSE3=ON -DFEELPP_TEST_SSE4_1=ON")
#   elseif(FEELPP_EXPLICIT_VECTORIZATION MATCHES SSE4_2)
#     set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEELPP_TEST_SSE2=ON -DFEELPP_TEST_SSE3=ON -DFEELPP_TEST_SSSE3=ON -DFEELPP_TEST_SSE4_1=ON -DFEELPP_TEST_SSE4_2=ON")
#   elseif(FEELPP_EXPLICIT_VECTORIZATION MATCHES Altivec)
#     set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEELPP_TEST_ALTIVEC=ON")
#   elseif(FEELPP_EXPLICIT_VECTORIZATION MATCHES novec)
#     set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEELPP_TEST_NO_EXPLICIT_VECTORIZATION=ON")
#   else(FEELPP_EXPLICIT_VECTORIZATION MATCHES SSE2)
#     message(FATAL_ERROR "Invalid value for FEELPP_EXPLICIT_VECTORIZATION (${FEELPP_EXPLICIT_VECTORIZATION}), must be: novec, SSE2, SSE3, Altivec")
#   endif(FEELPP_EXPLICIT_VECTORIZATION MATCHES SSE2)
# endif(DEFINED FEELPP_EXPLICIT_VECTORIZATION)

if(DEFINED FEELPP_CMAKE_ARGS)
  set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} ${FEELPP_CMAKE_ARGS}")
endif(DEFINED FEELPP_CMAKE_ARGS)

#The idea behind ctest launchers is that they wrap each compile or link step so
#the output can be saved and sent to CDash in the event of a warning or
#error. Rather than trying to grep through and analyze the full build output
#after thousands of compile and link calls, with this technique, ctest may
#simply capture the error output directly and pass it in its entirety to the
#dashboard. This helps immensely in figuring out some why some errors occur,
#without necessarily even having access to the client machine.
set(CTEST_USE_LAUNCHERS 1)


# raise the warning/error limit
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "33331")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS "33331")

# to get CTEST_PROJECT_SUBPROJECTS definition:
include("${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake")
# clear the binary directory and create an initial cache
#CTEST_EMPTY_BINARY_DIRECTORY (${CTEST_BINARY_DIRECTORY})
set(CTEST_INITIAL_CACHE "
CMAKE_CXX_COMPILER:STRING=${FEELPP_CXX}
FEELPP_ENABLE_ALL:BOOL=ON
")
# site
set(CTEST_SITE "${FEELPP_SITE}")
# build name
set(CTEST_BUILD_NAME "${FEELPP_BUILD_STRING}")
# should ctest wipe the binary tree before running
#SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

if(FEELPP_CXX AND NOT WIN32)
 set(CTEST_ENVIRONMENT "CXX=${FEELPP_CXX}")
endif(FEELPP_CXX AND NOT WIN32)
MESSAGE(WARNING "ctest_environment ${CTEST_ENVIRONMENT}")

MESSAGE(WARNING "FEELPP_MODE: ${FEELPP_MODE}")
ctest_start(${FEELPP_MODE})
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}")
ctest_submit(PARTS Update Notes)

message(WARNING "subprojects: ${CTEST_PROJECT_SUBPROJECTS}" )
foreach(subproject ${CTEST_PROJECT_SUBPROJECTS})
#foreach(subproject "opus")
  message(WARNING "testing subproject ${subproject}")
  set_property(GLOBAL PROPERTY SubProject ${subproject})
  set_property (GLOBAL PROPERTY Label ${subproject})
  ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND
    OPTIONS "-DCTEST_USE_LAUNCHERS=${CTEST_USE_LAUNCHERS};-DCMAKE_CXX_COMPILER:STRING=${FEELPP_CXX};-DFEELPP_ENABLE_ALL:BOOL=ON" )
  ctest_submit(PARTS Configure)
  message(WARNING "build target ${subproject}")
  #set(CTEST_BUILD_COMMAND "make ${FEELPP_MAKE_ARGS} -i ${subproject}")
  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND TARGET "${subproject}"  )
  # builds target ${CTEST_BUILD_TARGET}
  ctest_submit(PARTS Build)
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND INCLUDE_LABEL "${subproject}"  )
  # runs only tests that have a LABELS property matching "${subproject}"
  ctest_submit(PARTS Test)
endforeach()

