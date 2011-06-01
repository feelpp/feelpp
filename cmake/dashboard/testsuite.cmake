
####################################################################
#
# Usage:
#  - create a new folder, let's call it cdash
#  - in that folder, do:
#    ctest -S path/to/feel/cmake/dashboard/testsuite.cmake[,option1=value1[,option2=value2]]
#
# Options:
#  - FEEL_CXX: compiler, eg.: g++-4.2
#      default: default c++ compiler
#  - FEEL_SITE: eg, sd-25154, or the name of the contributor, etc.
#      default: hostname
#  - FEEL_BUILD_STRING: a string which identify the system/compiler. It should be formed like that:
#        <OS_name>-<OS_version>-<arch>-<compiler-version>
#      with:
#        <OS_name> = opensuse, debian, osx, windows, cygwin, freebsd, solaris, etc.
#        <OS_version> = 11.1, XP, vista, leopard, etc.
#        <arch> = i386, x86_64, ia64, powerpc, etc.
#        <compiler-version> = gcc-4.3.2, icc-11.0, MSVC-2008, etc.
#  - FEEL_EXPLICIT_VECTORIZATION: novec, SSE2, Altivec
#       default: SSE2 for x86_64 systems, novec otherwise
#       Its value is automatically appended to FEEL_BUILD_STRING
#  - FEEL_CMAKE_DIR: path to cmake executable
#  - FEEL_MODE: dashboard model, can be Experimental, Nightly, or Continuous
#      default: Nightly
#  - FEEL_WORK_DIR: directory used to download the source files and make the builds
#      default: folder which contains this script
#  - FEEL_CMAKE_ARGS: additional arguments passed to cmake
#  - FEEL_GENERATOR_TYPE: allows to overwrite the generator type
#      default: nmake (windows
#      See http://www.cmake.org/cmake/help/cmake2.6docs.html#section_Generators for a complete
#      list of supported generators.
#  - FEEL_NO_UPDATE: allows to submit dash boards from local repositories
#      This might be interesting in case you want to submit dashboards
#      including local changes.
#  - CTEST_SOURCE_DIRECTORY: path to feel's feel (use a new and empty folder, not the one you are working on)
#      default: <FEEL_WORK_DIR>/feel
#  - CTEST_BINARY_DIRECTORY: build directory
#      default: <FEEL_WORK_DIR>/nightly-<FEEL_CXX>
#
# Here is an example running several compilers on a linux system:
# #!/bin/bash
# ARCH=`uname -m`
# SITE=`hostname`
# VERSION=opensuse-11.1
# WORK_DIR=/home/prudhomm/scratch/cdash
# # get the last version of the script
# wget http://bitbucket.org/feel/feel/raw/tip/test/testsuite.cmake -o $WORK_DIR/testsuite.cmake
# COMMON="ctest -S $WORK_DIR/testsuite.cmake,FEEL_WORK_DIR=$WORK_DIR,FEEL_SITE=$SITE,FEEL_MODE=$1,FEEL_BUILD_STRING=$OS_VERSION-$ARCH"
# $COMMON-gcc-3.4.6,FEEL_CXX=g++-3.4
# $COMMON-gcc-4.0.1,FEEL_CXX=g++-4.0.1
# $COMMON-gcc-4.3.2,FEEL_CXX=g++-4.3,FEEL_EXPLICIT_VECTORIZATION=novec
# $COMMON-gcc-4.3.2,FEEL_CXX=g++-4.3,FEEL_EXPLICIT_VECTORIZATION=SSE2
# $COMMON-icc-11.0,FEEL_CXX=icpc
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

if(NOT FEEL_SITE)
  site_name(FEEL_SITE)
endif(NOT FEEL_SITE)

if(NOT FEEL_CMAKE_DIR)
  SET(FEEL_CMAKE_DIR "")
endif(NOT FEEL_CMAKE_DIR)

if (NOT FEEL_CXX)
  set(FEEL_CXX "g++")
endif(NOT FEEL_CXX)

if(NOT FEEL_BUILD_STRING)

  # let's try to find all information we need to make the build string ourself

  # OS
  build_name(FEEL_OS_VERSION)

  # arch
  set(FEEL_ARCH ${CMAKE_SYSTEM_PROCESSOR})
  if(WIN32)
    set(FEEL_ARCH $ENV{PROCESSOR_ARCHITECTURE})
  else(WIN32)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE FEEL_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif(WIN32)

  set(FEEL_BUILD_STRING ${FEEL_OS_VERSION}${FEEL_ARCH}-${FEEL_CXX})

endif(NOT FEEL_BUILD_STRING)

if(DEFINED FEEL_EXPLICIT_VECTORIZATION)
  set(FEEL_BUILD_STRING ${FEEL_BUILD_STRING}-${FEEL_EXPLICIT_VECTORIZATION})
endif(DEFINED FEEL_EXPLICIT_VECTORIZATION)

if(NOT FEEL_WORK_DIR)
  set(FEEL_WORK_DIR ${CTEST_SCRIPT_DIRECTORY})
endif(NOT FEEL_WORK_DIR)

if(NOT CTEST_SOURCE_DIRECTORY)
  SET (CTEST_SOURCE_DIRECTORY "${FEEL_WORK_DIR}/feel")
endif(NOT CTEST_SOURCE_DIRECTORY)

if(NOT CTEST_BINARY_DIRECTORY)
  SET (CTEST_BINARY_DIRECTORY "${FEEL_WORK_DIR}/nightly_${FEEL_CXX}")
endif(NOT CTEST_BINARY_DIRECTORY)

if(NOT FEEL_MODE)
  set(FEEL_MODE Nightly)
endif(NOT FEEL_MODE)

## mandatory variables (the default should be ok in most cases):

#if(NOT FEEL_NO_UPDATE)
SET (CTEST_SVN_COMMAND "svn")
SET (CTEST_SVN_CHECKOUT   "${CTEST_SVN_COMMAND} co svn://scm.forge.imag.fr/var/lib/gforge/chroot/scmrepos/svn/life/trunk/life/trunk ${CTEST_SOURCE_DIRECTORY}")
set (CTEST_UPDATE_COMMAND "${CTEST_SVN_COMMAND}")
  #SET(CTEST_BACKUP_AND_RESTORE TRUE) # the backup is SVN related ...
#endif(NOT FEEL_NO_UPDATE)

# which ctest command to use for running the dashboard
#SET (CTEST_COMMAND "${FEEL_CMAKE_DIR}ctest -D ${FEEL_MODE} --no-compress-output")
#if($ENV{FEEL_CTEST_ARGS})
#SET (CTEST_COMMAND "${CTEST_COMMAND} $ENV{FEEL_CTEST_ARGS}")
#endif($ENV{FEEL_CTEST_ARGS})
# what cmake command to use for configuring this dashboard
#SET (CTEST_CMAKE_COMMAND "${FEEL_CMAKE_DIR}cmake -DFEEL_LEAVE_TEST_IN_ALL_TARGET=ON")

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this make sure we get consistent outputs
SET($ENV{LC_MESSAGES} "en_EN")

# should ctest wipe the binary tree before running
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

# raise the warning/error limit
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "33331")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS "33331")

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
if(WIN32 AND NOT UNIX)
  #message(SEND_ERROR "win32")
  if(FEEL_GENERATOR_TYPE)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -G \"${FEEL_GENERATOR_TYPE}\"")
    SET (CTEST_INITIAL_CACHE "
      CMAKE_BUILD_TYPE:STRING=Release
      BUILDNAME:STRING=${FEEL_BUILD_STRING}
      SITE:STRING=${FEEL_SITE}
    ")
  else(FEEL_GENERATOR_TYPE)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -G \"NMake Makefiles\" -DCMAKE_MAKE_PROGRAM=nmake")
    SET (CTEST_INITIAL_CACHE "
      MAKECOMMAND:STRING=nmake /i
      CMAKE_MAKE_PROGRAM:FILEPATH=nmake
      CMAKE_GENERATOR:INTERNAL=NMake Makefiles
      CMAKE_BUILD_TYPE:STRING=Release
      BUILDNAME:STRING=${FEEL_BUILD_STRING}
      SITE:STRING=${FEEL_SITE}
    ")
  endif(FEEL_GENERATOR_TYPE)
else(WIN32 AND NOT UNIX)
  SET (CTEST_INITIAL_CACHE "
    BUILDNAME:STRING=${FEEL_BUILD_STRING}
    SITE:STRING=${FEEL_SITE}
  ")
endif(WIN32 AND NOT UNIX)

if (UNIX)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif(UNIX)

#set(CTEST_BUILD_COMMAND "${MAKE} ${OPTION_BUILD}")
set(CTEST_BUILD_COMMAND "make")
SET(CTEST_CMAKE_COMMAND "cmake" )
SET(CTEST_MAKE_COMMAND "${CMAKE_EXECUTABLE_NAME}" )
# set any extra environment variables to use during the execution of the script here:
# setting this variable on windows machines causes trouble ...

if(FEEL_CXX AND NOT WIN32)
  set(CTEST_ENVIRONMENT "CXX=${FEEL_CXX}")
endif(FEEL_CXX AND NOT WIN32)

if(DEFINED FEEL_EXPLICIT_VECTORIZATION)
  if(FEEL_EXPLICIT_VECTORIZATION MATCHES SSE2)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEEL_TEST_SSE2=ON")
  elseif(FEEL_EXPLICIT_VECTORIZATION MATCHES SSE3)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEEL_TEST_SSE2=ON -DFEEL_TEST_SSE3=ON")
  elseif(FEEL_EXPLICIT_VECTORIZATION MATCHES SSSE3)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEEL_TEST_SSE2=ON -DFEEL_TEST_SSE3=ON -DFEEL_TEST_SSSE3=ON")
  elseif(FEEL_EXPLICIT_VECTORIZATION MATCHES SSE4_1)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEEL_TEST_SSE2=ON -DFEEL_TEST_SSE3=ON -DFEEL_TEST_SSSE3=ON -DFEEL_TEST_SSE4_1=ON")
  elseif(FEEL_EXPLICIT_VECTORIZATION MATCHES SSE4_2)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEEL_TEST_SSE2=ON -DFEEL_TEST_SSE3=ON -DFEEL_TEST_SSSE3=ON -DFEEL_TEST_SSE4_1=ON -DFEEL_TEST_SSE4_2=ON")
  elseif(FEEL_EXPLICIT_VECTORIZATION MATCHES Altivec)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEEL_TEST_ALTIVEC=ON")
  elseif(FEEL_EXPLICIT_VECTORIZATION MATCHES novec)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DFEEL_TEST_NO_EXPLICIT_VECTORIZATION=ON")
  else(FEEL_EXPLICIT_VECTORIZATION MATCHES SSE2)
    message(FATAL_ERROR "Invalid value for FEEL_EXPLICIT_VECTORIZATION (${FEEL_EXPLICIT_VECTORIZATION}), must be: novec, SSE2, SSE3, Altivec")
  endif(FEEL_EXPLICIT_VECTORIZATION MATCHES SSE2)
endif(DEFINED FEEL_EXPLICIT_VECTORIZATION)

if(DEFINED FEEL_CMAKE_ARGS)
  set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} ${FEEL_CMAKE_ARGS}")
endif(DEFINED FEEL_CMAKE_ARGS)

#The idea behind ctest launchers is that they wrap each compile or link step so
#the output can be saved and sent to CDash in the event of a warning or
#error. Rather than trying to grep through and analyze the full build output
#after thousands of compile and link calls, with this technique, ctest may
#simply capture the error output directly and pass it in its entirety to the
#dashboard. This helps immensely in figuring out some why some errors occur,
#without necessarily even having access to the client machine.
set(CTEST_USE_LAUNCHERS 1)

# to get CTEST_PROJECT_SUBPROJECTS definition:
include("${CTEST_SOURCE_DIRECTORY}/../CTestConfig.cmake")

ctest_start(${FEEL_MODE})
#ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}")
#ctest_submit(PARTS Update Notes)
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" OPTIONS "-DCTEST_USE_LAUNCHERS=${CTEST_USE_LAUNCHERS}" APPEND)
#ctest_submit(PARTS Configure)

message(WARNING "subprojects: ${CTEST_PROJECT_SUBPROJECTS}" )
 foreach(subproject ${CTEST_PROJECT_SUBPROJECTS})
   message(WARNING "testing subproject ${subproject}")
   set_property(GLOBAL PROPERTY SubProject ${subproject})
   set_property (GLOBAL PROPERTY Label ${subproject})
   set(CTEST_BUILD_TARGET "${subproject}")
   message(WARNING "build target ${subproject}")
   ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND )
   # builds target ${CTEST_BUILD_TARGET}
# #  ctest_submit(PARTS Build)
#   ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" INCLUDE_LABEL "${subproject}" )
#   # runs only tests that have a LABELS property matching "${subproject}"
# #  ctest_submit(PARTS Test)
 endforeach()

